function ABCD_split_and_match(csvname, subj_list, race_hdr, site_hdr, fam_hdr, cfds_ls, bhvr_ls, niter, outdir, outname)

ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';

if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_list);
subj_hdr = 'subjectkey';

if(~exist('race_hdr', 'var') || isempty(race_hdr))
    race_hdr = 'race';
end

if(~exist('site_hdr', 'var') || isempty(site_hdr))
    site_hdr = 'site';
end
if(~exist('fam_hdr', 'var') || isempty(fam_hdr))
    fam_hdr = 'family_id';
end

if(~exist('cfds_ls', 'var') || isempty(cfds_hdr))
    cfds_ls = fullfile(ls_dir, 'confounds_list.txt');
end
cfds_nm = CBIG_text2cell(cfds_ls);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
bhvr_nm = CBIG_text2cell(bhvr_ls);

if(~exist('niter', 'var') || isempty(niter))
    niter = 100;
end

%% parse table
d = readtable(csvname);
row_idx = zeros(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects{s});
    if(all(tmp_idx==0))
        error('Subject %s does not exist in %s.', subjects{s}, csvname);
    else
        row_idx(s) = find(tmp_idx==1);
    end
end
d = d(row_idx,:);

race = d.(race_hdr);
site = d.(site_hdr);
for s = 1:nsub
    site{s} = str2num(site{s}(5:end));
end
site = cell2mat(site);

fam_id = d.(fam_hdr);

cfds = [];
n = strcmpi(cfds_nm, 'age');
if(any(n))
    age = d.(cfds_nm{n==1});
    cfds = [cfds age];
end

n = strcmpi(cfds_nm, 'sex');
if(any(n))
    sex = strcmp(d.(cfds_nm{n==1}), 'F');
    cfds = [cfds sex];
end

possib_cfds = {'FD', 'DVARS', 'ICV', 'peduc_avg'};
for c = 1:length(possib_cfds)
    n = strcmpi(cfds_nm, possib_cfds{c});
    if(any(n))
        col = d.(cfds_nm{n==1});
        cfds = [cfds col];
    end
end

cfds_zn = zscore(cfds, 0, 1);

bhvr = zeros(nsub, length(bhvr_nm));
for b = 1:length(bhvr_nm)
    bhvr(:,b) = d.(bhvr_nm{b});
end
bhvr_zn = zscore(bhvr, 0, 1);

%%
[selAA, selWA, sel_mAA, sel_mWA] = ABCD_match_WAtoAA_within_site(subjects, race, site, bhvr_zn, cfds_zn, niter);

for b = 1:length(bhvr_nm)
    curr_selAA = cat(2, selAA{b,:});
    [curr_selAA, curr_AAidx] = intersect(subjects, curr_selAA, 'stable');
    curr_AAsite = site(curr_AAidx);
    curr_AA_fam_id = fam_id(curr_AAidx);
    [fold_AA{b}, fold_sites{b}] = ABCD_split_folds(curr_selAA, curr_AAsite, curr_AA_fam_id, 9);
end


% fold_AA = ABCD_split_folds(subjects(race==2), site(race==2), fam_id(race==2));
% 
% selAA = cell(length(bhvr_nm), length(fold_AA));
% selWA = cell(length(bhvr_nm), length(fold_AA));
% sel_mAA = cell(length(bhvr_nm), 1);
% sel_mWA = cell(length(bhvr_nm), 1);
% for b = 1:length(bhvr_nm)
%     for f = 1:length(fold_AA)
%         [fold_AA{f}, AAidx_f] = intersect(subjects, fold_AA{f}, 'stable');
%         AA_sites_f = site(AAidx_f);
%         uq_AA_sites_f = unique(AA_sites_f);
%         
%         
%         for st = 1:length(uq_AA_sites_f)
%             possib_WA_idx = race==1 & site == uq_AA_sites_f(st);
%             
%             curr_AA = fold_AA{f}(AA_sites_f==uq_AA_sites_f(st));
%             curr_WA = subjects(possib_WA_idx)';
%             
%             if(~isempty(curr_WA))
%                 [~, curr_AAidx] = intersect(subjects, curr_AA, 'stable');
%                 [~, curr_WAidx] = intersect(subjects, curr_WA, 'stable');
%                 mAA = [cfds_zn(curr_AAidx, :) bhvr_zn(curr_AAidx, b)];
%                 mWA = [cfds_zn(curr_WAidx, :) bhvr_zn(curr_WAidx, b)];
%                 
%                 mAA_3d = reshape(mAA, size(mAA,1), 1, size(mAA,2));
%                 mWA_3d = reshape(mWA, 1, size(mWA,1), size(mWA,2));
%                 cost_mat = bsxfun(@minus, mAA_3d, mWA_3d);
%                 cost_mat = sum(abs(cost_mat),3);
%                 
%                 cost_last_iter = 1e4;
%                 curr_AA_newiter = curr_AA;
%                 mAA_newiter = mAA;
%                 for iter = 1:niter
%                     [asgn_WAidx_newiter, cost_new_iter] = munkres(cost_mat);
%                     while(any(asgn_WAidx_newiter==0))
%                         unmatch_idx = find(asgn_WAidx_newiter==0);
%                         curr_AA_newiter(unmatch_idx) = [];
%                         mAA_newiter(unmatch_idx,:) = [];
%                         cost_mat(unmatch_idx, :) = [];
%                         [asgn_WAidx_newiter, cost_new_iter] = munkres(cost_mat);
%                     end
%                     
%                     cost_new_iter = cost_new_iter / length(curr_AA_newiter);
%                     cost_idx = sub2ind(size(cost_mat), 1:1:length(curr_AA_newiter), asgn_WAidx_newiter);
%                     cost_currAA = cost_mat(cost_idx);
%                     [max_cost, maxidx] = max(cost_currAA);
%                     
%                     cost_diff = cost_last_iter - cost_new_iter;
%                     if(cost_diff / cost_new_iter <= 0.05 && max_cost<=3)
%                         break
%                     end
%                     
%                     curr_AA = curr_AA_newiter;
%                     asgn_WAidx = asgn_WAidx_newiter;
%                     mAA = mAA_newiter;
%                     cost_last_iter = cost_new_iter;
%                     
%                     curr_AA_newiter(maxidx) = [];
%                     cost_mat(maxidx,:) = [];
%                     mAA_newiter(maxidx,:) = [];
%                     if(length(curr_AA) == 1)
%                         if(cost_new_iter>3)
%                             curr_AA = [];
%                             mAA = [];
%                             asgn_WAidx = [];
%                             cost_currAA = [];
%                         end
%                         break
%                     end
%                 end
%                 
%                 disp(cost_currAA)
%                 selAA{b,f} = [selAA{b,f}; curr_AA];
%                 selWA{b,f} = [selWA{b,f}; curr_WA(asgn_WAidx)];
%                 sel_mAA{b} = [sel_mAA{b}; mAA];
%                 sel_mWA{b} = [sel_mWA{b}; mWA(asgn_WAidx,:)];
%             end
%         end
%     end
% end
% 
% fold_sz = zeros(length(fold_AA), 1);
% for f = 1:length(fold_AA)
%     [fold_AA{f}, AAidx_f] = intersect(subjects, fold_AA{f}, 'stable');
%     AA_sites_f = site(AAidx_f);
%     uq_AA_sites_f = unique(AA_sites_f);
%     
%     for st = 1:length(uq_AA_sites_f)
%     
%         fold_sz(f) = fold_sz(f) + length(find(site==uq_AA_sites_f(st)));
%     end
% end


end

