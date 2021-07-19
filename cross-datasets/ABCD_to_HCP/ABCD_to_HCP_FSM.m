function ABCD_to_HCP_FSM(ABCD_FC, HCP_FC, ABCD_subfold, ABCD_bhvr_nm, HCP_bhvr_nm, outdir)

% ABCD_to_HCP_FSM(ABCD_FC, HCP_FC, ABCD_subfold, ABCD_bhvr_nm, HCP_bhvr_nm, outdir)
%
% 

proj_dir = '/home/jingweil/storage/MyProject/fairAI';
if(~exist('HCP_FC', 'var') || isempty(HCP_FC))
    HCP_FC = fullfile(proj_dir, 'HCP_race', 'mat', 'RSFC_948.mat');
end
if(~exist('ABCD_FC', 'var') || isempty(ABCD_FC))
    ABCD_FC = fullfile(proj_dir, 'ABCD_race', 'mat', 'RSFC', 'pass_rs_pass_pheno_5351.mat');
end

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

%% prepare for vectorization of RSFC
temp = ones(nROI);
tril_ind = tril(temp, -1);
FC_hcp.corr_mat = reshape(FC_hcp.corr_mat, nROI^2, size(FC_hcp.corr_mat, 3));
FC_hcp.corr_mat = FC_hcp.corr_mat(tril_ind==1, :);
FC_abcd.corr_mat = reshape(FC_abcd.corr_mat, nROI^2, size(FC_abcd.corr_mat, 3));
FC_abcd.corr_mat = FC_abcd.corr_mat(tril_ind==1, :);

load(ABCD_subfold)
Nfolds = length(sub_fold);

%% compute for each ABCD fold
for f = 1:Nfolds
    % extract RSFC of ABCD training subjects
    ABCD_train = FC_abcd.corr_mat(:, sub_fold(f).fold_index==0);

    mkdir(fullfile(outdir, [ABCD_bhvr_nm '_' HCP_bhvr_nm], ['fold_' num2str(f)]))
    curr_out = fullfile(outdir, [ABCD_bhvr_nm '_' HCP_bhvr_nm], ['fold_' num2str(f)], ...
        'FSM_HCP_full_vs_ABCD.mat')
    if(~exist(curr_out, 'file'))
        FSM = CBIG_crossvalid_kernel_with_scale(ABCD_train, FC_hcp.corr_mat, [], [], 'corr');
        FSM = FSM( (size(ABCD_train,2)+1):end, 1:size(ABCD_train,2) );
        save(curr_out, 'FSM')
    end
end
    
end