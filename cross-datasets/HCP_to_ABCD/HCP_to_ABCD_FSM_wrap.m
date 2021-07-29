function HCP_to_ABCD_FSM_wrap(outdir, HCP_FC, ABCD_FC, HCP_bhvr_ls, ABCD_bhvr_ls, ...
    ABCD_full_bhvr_ls, ABCD_subj_ls, HCP_split_dir, max_HCP_seed, ABCD_selAAWA)

% HCP_to_ABCD_FSM_wrap(outdir, HCP_FC, ABCD_FC, HCP_bhvr_ls, ABCD_bhvr_ls, ABCD_full_bhvr_ls, ABCD_subj_ls)
%
% 

proj_dir = '/home/jingweil/storage/MyProject/fairAI';
if(~exist('HCP_FC', 'var') || isempty(HCP_FC))
    HCP_FC = fullfile(proj_dir, 'HCP_race', 'mat', 'RSFC_948.mat');
end
if(~exist('ABCD_FC', 'var') || isempty(ABCD_FC))
    ABCD_FC = fullfile(proj_dir, 'ABCD_race', 'mat', 'RSFC', 'pass_rs_pass_pheno_5351_z.mat');
end
if(~exist('HCP_bhvr_ls', 'var') || isempty(HCP_bhvr_ls))
    HCP_bhvr_ls = fullfile(proj_dir, 'cross_ABCD_HCP', 'scripts', 'lists', 'HCP_behavior.txt');
end
[HCP_bhvrs, nbhvr] = CBIG_text2cell(HCP_bhvr_ls);
if(~exist('ABCD_bhvr_ls', 'var') || isempty(ABCD_bhvr_ls))
    ABCD_bhvr_ls = fullfile(proj_dir, 'cross_ABCD_HCP', 'scripts', 'lists', 'ABCD_behavior.txt');
end
ABCD_bhvrs = CBIG_text2cell(ABCD_bhvr_ls);
if(~exist('ABCD_full_bhvr_ls', 'var') || isempty(ABCD_full_bhvr_ls))
    ABCD_full_bhvr_ls = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'behavior_list.txt');
end
if(~exist('ABCD_subj_ls', 'var') || isempty(ABCD_subj_ls))
    ABCD_subj_ls = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt');
end

parfor b = 1:nbhvr
    HCP_bhvr_nm = HCP_bhvrs{b};
    ABCD_bhvr_nm = ABCD_bhvrs{b};
    fprintf('%s -- %s\n', HCP_bhvr_nm, ABCD_bhvr_nm)

    for seed = 1:max_HCP_seed
        HCP_subfold = fullfile(HCP_split_dir, ['split_seed' num2str(seed)], ...
            ['no_relative_10_fold_sub_list_' HCP_bhvr_nm '.mat']);
        if(~exist(HCP_subfold, 'file')); continue; end
        fprintf('\tSeed: %d\n', seed)
        HCP_to_ABCD_FSM(HCP_FC, ABCD_FC, HCP_subfold, ABCD_selAAWA, HCP_bhvr_nm, ABCD_bhvr_nm, ...
            fullfile(outdir, ['randseed_' num2str(seed)]), ABCD_subj_ls, ABCD_full_bhvr_ls);
    end
end

end
