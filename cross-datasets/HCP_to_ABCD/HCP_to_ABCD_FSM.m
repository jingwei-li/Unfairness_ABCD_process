function HCP_to_ABCD_FSM(HCP_FC, ABCD_FC, HCP_subfold, ABCD_selAAWA, HCP_bhvr_nm, ABCD_bhvr_nm, ...
    outdir, ABCD_subj_ls, ABCD_bhvr_ls)

% HCP_to_ABCD_FSM(HCP_FC, ABCD_FC, HCP_subfold, ABCD_selAAWA, HCP_bhvr_nm, ABCD_bhvr_nm, ...
%    outdir, ABCD_subj_ls, HCP_bhvr_ls, ABCD_bhvr_ls)
%
% Compute the functional similarity matrix between the training subjects in the HCP datasets
% and the test subjects in the ABCD datasets.
%
% Inputs:
% - HCP_FC (optional)
%   The full path of the RSFC file among all HCP subjects.
%   Default: /home/jingweil/storage/MyProject/fairAI/HCP_race/mat/RSFC_948.mat
% - ABCD_FC (optional)
%   The full path of the RSFC file among all ABCD subjects.
%   Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/mat/RSFC/pass_rs_pass_pheno_5351.mat
% - HCP_subfold
%   The full path of the split sub-folds file of HCP dataset for the current behavioral measure.
% - ABCD_selAAWA
%   The full path of the .mat file which contains the selected AA and WA.
% - HCP_bhvr_nm
%   The behavioral name used in the HCP dataset of the current NIH measure under investigation.
% - ABCD_bhvr_nm
%   The behavioral name used in the ABCD dataset of the current NIH measure under investigation.
% - outdir
% - ABCD_subj_ls (optional)
%   The full path of the subject list with all ABCD participants involved in this project.
%   It is used to index the selected AA and WA among all subjects, to extract the corresponding RSFC.
%   Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt
% - ABCD_bhvr_ls (optional)
%   The full path of all ABCD behavioral names. It is used to index the current behavioral measure
%   from all measures, so that the selected AA and WA (specifically for this behavior) can be identified
%   from `ABCD_selAAWA`.
%   Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt
%
% Author: Jingwei Li

proj_dir = '/home/jingweil/storage/MyProject/fairAI';
if(isempty(HCP_FC))
    HCP_FC = fullfile(proj_dir, 'HCP_race', 'mat', 'RSFC_948.mat');
end
if(~exist('ABCD_FC', 'var') || isempty(ABCD_FC))
    ABCD_FC = fullfile(proj_dir, 'ABCD_race', 'mat', 'RSFC', 'pass_rs_pass_pheno_5351.mat');
end
if(~exist('ABCD_subj_ls', 'var') || isempty(ABCD_subj_ls))
    ABCD_subj_ls = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt');
end
[all_ABCD_subj, nABCD] = CBIG_text2cell(ABCD_subj_ls);
if(~exist('ABCD_bhvr_ls', 'var') || isempty(ABCD_bhvr_ls))
    ABCD_bhvr_ls = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'behavior_list.txt');
end
all_ABCD_bhvr = CBIG_text2cell(ABCD_bhvr_ls);

%% load full RSFC matrices of both datasets
if(ischar(HCP_FC))
    FC_hcp = load(HCP_FC);
else
    FC_hcp.corr_mat = HCP_FC;
    clear HCP_FC
end
if(ischar(ABCD_FC))
    FC_abcd = load(ABCD_FC);
else
    FC_abcd.corr_mat = ABCD_FC;
end
nROI = size(FC_hcp.corr_mat, 1);

%% load data split of HCP subjects for the current behavior
load(HCP_subfold)
Nfolds = length(sub_fold);

%% extract selected AA and WA in the ABCD dataset for the current behavior
ABCD_sel = load(ABCD_selAAWA);
[~, bidx] = intersect(all_ABCD_bhvr, ABCD_bhvr_nm, 'stable');
ABCD_selAA = cat(2, ABCD_sel.selAA{bidx, :});
ABCD_selWA = cat(1, ABCD_sel.selWA{bidx, :});
[ABCD_selAA, AA_idx] = intersect(all_ABCD_subj, ABCD_selAA, 'stable');
[ABCD_selWA, WA_idx] = intersect(all_ABCD_subj, ABCD_selWA, 'stable');

%% prepare for vectorization of RSFC
temp = ones(nROI);
tril_ind = tril(temp, -1);
FC_abcd.corr_mat = reshape(FC_abcd.corr_mat, nROI^2, nABCD);
FC_abcd.corr_mat = FC_abcd.corr_mat(tril_ind==1, :);

%% compute for each HCP fold
for f = 1:Nfolds
    % extract RSFC of HCP training subjects
    HCP_train = FC_hcp.corr_mat(:,:, sub_fold(f).fold_index==0);
    HCP_train = reshape(HCP_train, nROI^2, size(HCP_train, 3));
    HCP_train = HCP_train(tril_ind==1, :);

    % extract RSFC of ABCD AA & WA separately
    ABCD_testAA = FC_abcd.corr_mat(:, AA_idx);
    ABCD_testWA = FC_abcd.corr_mat(:, WA_idx);

    % compute FSM
    mkdir(fullfile(outdir, [HCP_bhvr_nm '_' ABCD_bhvr_nm], ['fold_' num2str(f)]))
    curr_out = fullfile(outdir, [HCP_bhvr_nm '_' ABCD_bhvr_nm], ...
        ['fold_' num2str(f)], 'FSM_ABCD_full_vs_HCP.mat');
    if(~exist(curr_out, 'file'))
        FSM = CBIG_crossvalid_kernel_with_scale(HCP_train, FC_abcd.corr_mat, [], [], 'corr');
        FSM = FSM( (size(HCP_train,2)+1):end, 1:size(HCP_train,2) );
        save(curr_out, 'FSM')
    end

    curr_out = fullfile(outdir, [HCP_bhvr_nm '_' ABCD_bhvr_nm], ...
        ['fold_' num2str(f)], 'FSM_ABCD_AA_vs_HCP.mat');
    if(~exist(curr_out, 'file'))
        FSM = CBIG_crossvalid_kernel_with_scale(HCP_train, ABCD_testAA, [], [], 'corr');
        FSM = FSM( (size(HCP_train,2)+1):end, 1:size(HCP_train,2) );
        save(curr_out, 'FSM')
    end

    curr_out = fullfile(outdir, [HCP_bhvr_nm '_' ABCD_bhvr_nm], ...
        ['fold_' num2str(f)], 'FSM_ABCD_WA_vs_HCP.mat');
    if(~exist(curr_out, 'file'))
        FSM = CBIG_crossvalid_kernel_with_scale(HCP_train, ABCD_testWA, [], [], 'corr');
        FSM = FSM( (size(HCP_train,2)+1):end, 1:size(HCP_train,2) );
        save(curr_out, 'FSM')
    end
end

end
