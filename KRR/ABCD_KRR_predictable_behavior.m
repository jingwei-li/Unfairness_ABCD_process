function ABCD_KRR_predictable_behavior(bhvr_ls, colloq_ls, model_dir, split_dir, split_fstem, metric, outmat)

% ABCD_KRR_predictable_behavior(bhvr_ls, model_dir, split_dir, split_fstem, metric, outmat)
% 
% Perform permutation test (shullfing predicted scores) to obtain behaviors
% with statistically significant accuracy.
%
% Example:
% ABCD_KRR_predictable_behavior([], [], ...
%    '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/models/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y', ...
%    '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/mat/matchANDsplit/20200719', ...
%    '_pass_rs_pass_pheno', 'predictive_COD', ...
%    '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/mat/predictability/pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y.mat')

ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
colloq_nm = CBIG_text2cell(colloq_ls);

rng('default'); rng(1);
nperm = 1000;
alpha = 0.05;

null_stats = cell(nbhvr, 1);
p_perm = nan(nbhvr, 1);
for b = 1:nbhvr
    fprintf('#%d behavior: %s ...\n', b, bhvr_nm{b})
    fold_file = fullfile(split_dir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']);
    load(fold_file)
    krr_file = fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']);
    krr = load(krr_file);
    
    %% load original scores and predicted scores, and compute null distributions
    null_stats{b} = zeros(length(sub_fold), nperm);
    for f = 1:length(sub_fold)
        krry = load(fullfile(model_dir, 'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']));
        test_idx = sub_fold(f).fold_index==0;
        y_test = krry.y_resid(test_idx);
        y_train = krry.y_resid(~test_idx);
        y_pred = krr.y_predict_concat(test_idx);
        
        % compute null distributions
        for n = 1:nperm
            null_idx = datasample(1:length(y_test), length(y_test), 'replace', false);
            null_y_pred = y_pred(null_idx);
            
            null_stats{b}(f,n) = CBIG_compute_prediction_acc_and_loss(null_y_pred, y_test, ...
                metric, y_train);
        end
    end
    
    %% compute p values
    avg_stats = mean(krr.optimal_stats.(metric), 1);
    avg_null_stats = mean(null_stats{b}, 1);
    p_perm(b) = length(find(avg_null_stats > avg_stats)) ./ nperm;
    
    clear sub_fold krr
end

H_perm_FDR = FDR(p_perm, alpha);
sig_perm_idx = sort(H_perm_FDR);
sig_perm_bhvr = bhvr_nm(sig_perm_idx);
sig_perm_colloq = colloq_nm(sig_perm_idx);

%% output
outdir = fileparts(outmat);
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(outmat, 'null_stats', 'alpha', 'p_perm', 'H_perm_FDR', 'sig_perm_idx', 'sig_perm_bhvr', ...
    'sig_perm_colloq')

mkdir(fullfile(model_dir, 'lists'))
outtxt = fullfile(model_dir, 'lists', ['predictable_behaviors_' metric '.txt']);
% CBIG_cell2text(sig_perm_bhvr, outtxt);
% outtxt = fullfile(model_dir, 'lists', ['predictable_colloquial_' metric '.txt']);
% CBIG_cell2text(sig_perm_colloq, outtxt);
    
end