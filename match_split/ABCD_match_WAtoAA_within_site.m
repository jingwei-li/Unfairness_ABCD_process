function [selAA, selWA, sel_mAA, sel_mWA] = ABCD_match_WAtoAA_within_site( ...
    subjects, race, site, bhvr_zn, cfds_zn, niter, cost_ceil)

% [selAA, selWA, sel_mAA, sel_mWA] = ABCD_match_WAtoAA_within_site( ...
%     subjects, race, site, bhvr_zn, cfds_zn, niter, cost_ceil)
%
% For each single site, match African Americans (AA) with white Americans (WA) 
% using Hungarian matching.
%
% Inputs:
% - subjects
%   A cell array containing all subject IDs.
%
% - race
%   A column vector with the ethnicity/race label of each subject.
%
% - site
%   A column vector with the site ID of each subject.
%
% - bhvr_zn
%   A #subjects x #behaviors matrix, z-score normalized behavioral scores.
%
% - cfds_zn
%   A #subjects x #confounds matrix, z-score normalized confounding variables.
% 
% - niter
%   Maximal number of iterations used to exclude AA with highest matching cost.
% 
% - cost_ceil
%   The upper bound for matching cost. `niter` and `cost_ceil` work together to 
%   stop the iterations of subjects exclusion.
%
% Outputs:
% - selAA
%   A #behaviors x #sites cell. Each entry is the selected AA IDs which can be 
%   matched for a certain behavioral measure within a certain site.
%
% - selWA 
%   A #behaviors x #sites cell. Each entry is the selected WA IDs which can be 
%   matched for a certain behavioral measure within a certain site.
%
% - sel_mAA
%   A cell array with a length: number of behavioral measures. Each entry is a 
%   matrix (dimension: #confounds+1 x #AA) containing the confounds and the scores of
%   a certain behavioral measure of all matched AA concatenated across sites.
%
% - sel_mWA
%   A cell array with a length: number of behavioral measures. Each entry is a 
%   matrix (dimension: #confounds+1 x #WA) containing the confounds and the scores of
%   a certain behavioral measure of all matched WA concatenated across sites.
%
% Author: Jingwei Li

WArace = 1;
AArace = 2;
if(~exist('cost_ceil', 'var') || isempty(cost_ceil))
    cost_ceil = 2.45;
end

AA = subjects(race == AArace);
AA_sites = site(race == AArace);
uq_AA_sites = unique(AA_sites);

selAA = cell(size(bhvr_zn,2), length(uq_AA_sites));
selWA = cell(size(bhvr_zn,2), length(uq_AA_sites));
sel_mAA = cell(size(bhvr_zn,2), 1);
sel_mWA = cell(size(bhvr_zn,2), 1);
for b = 1:size(bhvr_zn,2)
    fprintf('Behavior #%d ...\n', b);
    for st = 1:length(uq_AA_sites)
        possib_WA_idx = race == WArace & site == uq_AA_sites(st);
        
        curr_AA = AA(AA_sites==uq_AA_sites(st));
        curr_WA = subjects(possib_WA_idx)';
        
        if(~isempty(curr_WA))
            [~, curr_AAidx] = intersect(subjects, curr_AA, 'stable');
            [~, curr_WAidx] = intersect(subjects, curr_WA, 'stable');

            mAA = bhvr_zn(curr_AAidx, b); mWA = bhvr_zn(curr_WAidx, b);
            if(~isempty(cfds_zn))
                mAA = [cfds_zn(curr_AAidx, :) mAA];
                mWA = [cfds_zn(curr_WAidx, :) mWA];
            end
            
            mAA_3d = reshape(mAA, size(mAA,1), 1, size(mAA,2));
            mWA_3d = reshape(mWA, 1, size(mWA,1), size(mWA,2));
            cost_mat = bsxfun(@minus, mAA_3d, mWA_3d);
            cost_mat = sum(abs(cost_mat),3);
            
            cost_last_iter = 1e4;
            curr_AA_newiter = curr_AA;
            mAA_newiter = mAA;
            for iter = 1:niter
                [asgn_WAidx_newiter, cost_new_iter] = munkres(cost_mat);
                while(any(asgn_WAidx_newiter==0))
                    unmatch_idx = find(asgn_WAidx_newiter==0);
                    curr_AA_newiter(unmatch_idx) = [];
                    mAA_newiter(unmatch_idx,:) = [];
                    cost_mat(unmatch_idx, :) = [];
                    [asgn_WAidx_newiter, cost_new_iter] = munkres(cost_mat);
                    curr_AA = curr_AA_newiter;
                end
                
                cost_new_iter = cost_new_iter / length(curr_AA_newiter);
                cost_idx = sub2ind(size(cost_mat), 1:1:length(curr_AA_newiter), asgn_WAidx_newiter);
                cost_currAA = cost_mat(cost_idx);
                [max_cost, maxidx] = max(cost_currAA);
                
                if(iter==1)
                    asgn_WAidx = asgn_WAidx_newiter;
                end
                
                cost_diff = cost_last_iter - cost_new_iter;
                if(cost_new_iter<1e-10 || (cost_diff / cost_new_iter <= 0.05 && max_cost<=cost_ceil))
                    break
                end
                
                curr_AA = curr_AA_newiter;
                asgn_WAidx = asgn_WAidx_newiter;
                mAA = mAA_newiter;
                cost_last_iter = cost_new_iter;
                
                curr_AA_newiter(maxidx) = [];
                cost_mat(maxidx,:) = [];
                mAA_newiter(maxidx,:) = [];
                if(length(curr_AA) == 1)
                    if(cost_new_iter>cost_ceil)
                        curr_AA = [];
                        mAA = [];
                        asgn_WAidx = [];
                        cost_currAA = [];
                    end
                    break
                end
            end
            
            %disp(cost_currAA)
            selAA{b,st} = [selAA{b,st}; curr_AA];
            selWA{b,st} = [selWA{b,st}; curr_WA(asgn_WAidx)];
            sel_mAA{b} = [sel_mAA{b}; mAA];
            sel_mWA{b} = [sel_mWA{b}; mWA(asgn_WAidx,:)];
        end
    end
end

end

