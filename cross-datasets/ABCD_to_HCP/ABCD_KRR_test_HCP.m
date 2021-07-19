function ABCD_KRR_test_HCP(outdir, ABCD_KRR_dir, ABCD_split_dir, crs_dt_FSM_dir, ABCD_bhvr_nm, HCP_bhvr_nm, ...
    HCP_KRR_dir, max_HCP_seed, HCP_splitAAdir, HCP_splitWAdir, usable_seeds_lsdir, HCP_subj_ls)

% ABCD_KRR_test_HCP(outdir, ABCD_KRR_dir, crs_dt_FSM_dir, ABCD_bhvr_nm, HCP_bhvr_nm, ...
%     HCP_KRR_dir, max_HCP_seed, HCP_splitAAdir, HCP_splitWAdir, usable_seeds_lsdir)
%
% 

proj_dir = '/home/jingweil/storage/MyProject/fairAI/';
HCP_58_ls = fullfile(proj_dir, 'HCP_race', 'scripts', 'lists', ...
    'Cognitive_Personality_Task_Social_Emotion_58.txt');
HCP_58 = CBIG_text2cell(HCP_58_ls);
b_idx = strcmp(HCP_58, HCP_bhvr_nm);

[HCP_subj] = CBIG_text2cell(HCP_subj_ls);
HCP_y = load(fullfile(HCP_KRR_dir, ['y_' HCP_bhvr_nm '.mat']));

bin_flag = 0;
with_bias = 1;
metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};


load(fullfile(ABCD_split_dir, ['sub_fold_pass_rs_pass_pheno_' ABCD_bhvr_nm '.mat']));
Nfolds = length(sub_fold);
opt_params = load(fullfile(ABCD_KRR_dir, ABCD_bhvr_nm, ['final_result_' ABCD_bhvr_nm '.mat']));
for f = 1:Nfolds
    ABCD_y = load(fullfile(ABCD_KRR_dir, ABCD_bhvr_nm, 'y', ['fold_' num2str(f)], ['y_regress_' ABCD_bhvr_nm '.mat']));
    y_train{f} = ABCD_y.y_resid(sub_fold(f).fold_index==0);
    FSM_train = load(fullfile(ABCD_KRR_dir, ABCD_bhvr_nm, 'FSM_innerloop', ['fold_' num2str(f)], 'FSM_corr.mat'));
    FSM_hcp = load(fullfile(crs_dt_FSM_dir, [ABCD_bhvr_nm '_' HCP_bhvr_nm], ...
        ['fold_' num2str(f)], 'FSM_HCP_full_vs_ABCD.mat'));
    [y_p{f}, y_t{f}, acc(f), pred_stats{f}] = CBIG_KRR_test_cv( bin_flag, ...
        FSM_train.FSM, FSM_hcp.FSM, y_train{f}, HCP_y.y, HCP_y.y, ...
        with_bias, opt_params.optimal_lambda(f), opt_params.optimal_threshold(f), metrics );

    c = 0;
    for seed = 1:max_HCP_seed
        usable_txt = fullfile(usable_seeds_lsdir, ['usable_behaviors_seed' num2str(seed) '.txt']);
        if(exist(usable_txt, 'file'))
            mtch_bhvr = CBIG_text2cell(usable_txt);
            idx = strcmp(mtch_bhvr, HCP_bhvr_nm);
            if(any(idx==1))
                c = c+1;
                load(fullfile(HCP_splitAAdir, ['split_seed' num2str(seed) '.mat']))
                if(exist('Afr_fold', 'var'))
                    AA_fold = Afr_fold; 
                    clear Afr_fold
                end
                load(fullfile(HCP_splitWAdir, ['split_seed' num2str(seed), '.mat']))

                [~, AA_idx] = intersect(HCP_subj, cat(1, AA_fold.sub_perfold{:}), 'stable');
                [~, WA_idx] = intersect(HCP_subj, cat(1, best_assign{b_idx}{:}), 'stable');
                y_p_AA{f,c} = y_p{f}{1}(AA_idx);   y_t_AA{f,c} = y_t{f}{1}(AA_idx);
                y_p_WA{f,c} = y_p{f}{1}(WA_idx);   y_t_WA{f,c} = y_t{f}{1}(WA_idx);
                acc_AA(f,c) = CBIG_corr(y_p_AA{f,c}, y_t_AA{f,c});
                acc_WA(f,c) = CBIG_corr(y_p_WA{f,c}, y_t_WA{f,c});
                for k = 1:length(metrics)
                    [pred_stats_AA{f,c}(k,1),~] = CBIG_compute_prediction_acc_and_loss(...
                        y_p_AA{f,c}, y_t_AA{f,c}, metrics{k}, y_train{f});
                    [pred_stats_WA{f,c}(k,1),~] = CBIG_compute_prediction_acc_and_loss(...
                        y_p_WA{f,c}, y_t_WA{f,c}, metrics{k}, y_train{f});
                end
            end
        end
    end
end
    
outmat = fullfile(outdir, ABCD_bhvr_nm, 'prediction.mat');
mkdir(fullfile(outdir, ABCD_bhvr_nm))
save(outmat, 'y_p', 'y_t', 'y_train', 'acc', 'pred_stats', ...
    'y_p_AA', 'y_t_AA', 'acc_AA', 'pred_stats_AA', ...
    'y_p_WA', 'y_t_WA', 'acc_WA', 'pred_stats_WA')
end