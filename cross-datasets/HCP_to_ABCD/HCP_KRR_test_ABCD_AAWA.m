function HCP_KRR_test_ABCD_AAWA(outdir, HCP_KRR_dir, max_HCP_seed, ABCD_selAAWA, HCP_bhvr_nm, ABCD_bhvr_nm, ...
    crs_dt_FSM_dir, ABCD_csv, ABCD_subj_ls, ABCD_bhvr_ls, HCP_splitAAdir, HCP_splitWAdir)

% HCP_KRR_test_ABCD_AAWA(outdir, HCP_KRR_dir, max_HCP_seed, ABCD_selAAWA, HCP_bhvr_nm, ABCD_bhvr_nm, ...
%     crs_dt_FSM_dir, ABCD_csv, ABCD_subj_ls, HCP_splitAAdir, HCP_splitWAdir)
%
% 

proj_dir = '/home/jingweil/storage/MyProject/fairAI/';
if(~exist('ABCD_csv', 'var') || isempty(ABCD_csv))
    ABCD_csv = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'phenotypes_pass_rs.txt');
end
if(~exist('ABCD_subj_ls', 'var') || isempty(ABCD_subj_ls))
    ABCD_subj_ls = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt');
end
if(~exist('ABCD_bhvr_ls', 'var') || isempty(ABCD_bhvr_ls))
    ABCD_bhvr_ls = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'behavior_list.txt');
end
HCP_58_ls = fullfile(proj_dir, 'HCP_race', 'scripts', 'lists', ...
    'Cognitive_Personality_Task_Social_Emotion_58.txt');
HCP_58 = CBIG_text2cell(HCP_58_ls);
b_idx = strcmp(HCP_58, HCP_bhvr_nm);

ABCD_subj = CBIG_text2cell(ABCD_subj_ls);
d_abcd = readtable(ABCD_csv);
[~, ~, idx] = intersect(ABCD_subj, d_abcd.subjectkey, 'stable');
ABCD_y = d_abcd.(ABCD_bhvr_nm)(idx);

all_ABCD_bhvr = CBIG_text2cell(ABCD_bhvr_ls);
[~, bidx] = intersect(all_ABCD_bhvr, ABCD_bhvr_nm, 'stable');
ABCD_sel = load(ABCD_selAAWA);
ABCD_selAA = cat(2, ABCD_sel.selAA{bidx, :});
ABCD_selWA = cat(1, ABCD_sel.selWA{bidx, :});
[ABCD_selAA, AA_idx] = intersect(ABCD_subj, ABCD_selAA, 'stable');
[ABCD_selWA, WA_idx] = intersect(ABCD_subj, ABCD_selWA, 'stable');

y_testAA = ABCD_y(AA_idx);
y_testWA = ABCD_y(WA_idx);

bin_flag = 0;
with_bias = 1;
metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};

for seed = 1:max_HCP_seed
    curr_model_dir = fullfile(HCP_KRR_dir, ['randseed_' num2str(seed)], HCP_bhvr_nm);
    if(~exist(curr_model_dir, 'dir')); continue; end

    load(fullfile(HCP_splitAAdir, ['split_seed' num2str(seed) '.mat']))
    %%%%%% for author's own debug purpose
    if(~exist('AA_fold', 'var'))
        AA_fold = Afr_fold; clear Afr_fold
    end
    load(fullfile(HCP_splitWAdir, ['split_seed' num2str(seed) '.mat']))

    curr_outdir = fullfile(outdir, ['randseed_' num2str(seed)], HCP_bhvr_nm);
    mkdir(curr_outdir)
    opt_params = load(fullfile(curr_model_dir, ['final_result_' HCP_bhvr_nm '.mat']));

    load(fullfile(curr_model_dir, ['no_relative_10_fold_sub_list_' HCP_bhvr_nm '.mat']))
    Nfolds = length(sub_fold);
    for f = 1:Nfolds
        HCP_y = load(fullfile(curr_model_dir, 'y', ['fold_' num2str(f)], ['y_regress_' HCP_bhvr_nm '.mat']));
        y_train{f} = HCP_y.y_resid(sub_fold(f).fold_index==0);
        AA_train = cat(1, AA_fold.sub_perfold{setdiff(1:Nfolds, f)});
        [~, idx_AAtrain] = intersect(sub_fold(1).all_subjects, AA_train, 'stable');
        y_trainAA{f} = HCP_y.y_resid(idx_AAtrain);
        WA_train = cat(1, best_assign{b_idx}{setdiff(1:Nfolds, f)});
        [~, idx_WAtrain] = intersect(sub_fold(1).all_subjects, WA_train, 'stable');
        y_trainWA{f} = HCP_y.y_resid(idx_WAtrain);

        FSM_train = load(fullfile(curr_model_dir, 'FSM_innerloop', ['fold_' num2str(f)], 'FSM_corr.mat'));
        FSM_testAA = load(fullfile(crs_dt_FSM_dir, ['randseed_' num2str(seed)], [HCP_bhvr_nm '_' ABCD_bhvr_nm], ...
            ['fold_' num2str(f)], 'FSM_ABCD_AA_vs_HCP.mat'));
        FSM_testWA = load(fullfile(crs_dt_FSM_dir, ['randseed_' num2str(seed)], [HCP_bhvr_nm '_' ABCD_bhvr_nm], ...
            ['fold_' num2str(f)], 'FSM_ABCD_WA_vs_HCP.mat'));
        FSM_abcd = load(fullfile(crs_dt_FSM_dir, ['randseed_' num2str(seed)], [HCP_bhvr_nm '_' ABCD_bhvr_nm], ...
            ['fold_' num2str(f)], 'FSM_ABCD_full_vs_HCP.mat'));

        [y_p{f}, y_t{f}, acc(f), pred_stats{f}] = CBIG_KRR_test_cv( bin_flag, ...
            FSM_train.FSM, FSM_abcd.FSM, y_train{f}, ABCD_y, ABCD_y, ...
            with_bias, opt_params.optimal_lambda(f), opt_params.optimal_threshold(f), metrics );

        [y_p_AA{f}, y_t_AA{f}, acc_AA(f), pred_stats_AA{f}] = CBIG_KRR_test_cv( bin_flag, ...
            FSM_train.FSM, FSM_testAA.FSM, y_train{f}, y_testAA, y_testAA, ...
            with_bias, opt_params.optimal_lambda(f), opt_params.optimal_threshold(f), metrics );
        
        [y_p_WA{f}, y_t_WA{f}, acc_WA(f), pred_stats_WA{f}] = CBIG_KRR_test_cv( bin_flag, ...
            FSM_train.FSM, FSM_testWA.FSM, y_train{f}, y_testWA, y_testWA, ...
            with_bias, opt_params.optimal_lambda(f), opt_params.optimal_threshold(f), metrics );
    end
    save(fullfile(curr_outdir, 'prediction_AAWA.mat'), 'y_p', 'y_t', 'acc', 'pred_stats', ...
        'y_p_AA', 'y_p_WA', 'y_t_AA', 'y_t_WA', 'acc_AA', 'acc_WA', 'pred_stats_AA', ...
        'pred_stats_WA', 'y_train', 'y_trainAA', 'y_trainWA')
end