function HCP_KRR_test_ABCD_predictable_behavior(outmat, train_dir, test_dir, testFSM_dir, max_HCP_seed, Nperm, test_metric, ...
    HCP_bhvr_ls, HCP_colloq_ls, ABCD_bhvr_ls, ABCD_colloq_ls, ABCD_csv, ABCD_subj_ls)

% HCP_KRR_test_ABCD_predictable_behavior(input)
%
% 

%% default filenames
proj_dir = '/home/jingweil/storage/MyProject/fairAI';
if(~exist('HCP_bhvr_ls', 'var') || isempty(HCP_bhvr_ls))
    HCP_bhvr_ls = fullfile(proj_dir, 'cross_ABCD_HCP', 'scripts', 'lists', 'HCP_behavior.txt');
end
if(~exist('HCP_colloq_ls', 'var') || isempty(HCP_colloq_ls))
    HCP_colloq_ls = fullfile(proj_dir, 'cross_ABCD_HCP', 'scripts', 'lists', 'HCP_colloquial.txt');
end
if(~exist('ABCD_bhvr_ls', 'var') || isempty(ABCD_bhvr_ls))
    ABCD_bhvr_ls = fullfile(proj_dir, 'cross_ABCD_HCP', 'scripts', 'lists', 'ABCD_behavior.txt');
end
if(~exist('ABCD_colloq_ls', 'var') || isempty(ABCD_colloq_ls))
    ABCD_colloq_ls = fullfile(proj_dir, 'cross_ABCD_HCP', 'scripts', 'lists', 'ABCD_colloquial.txt');
end
if(~exist('ABCD_csv', 'var') || isempty(ABCD_csv))
    ABCD_csv = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'phenotypes_pass_rs.txt');
end
if(~exist('ABCD_subj_ls', 'var') || isempty(ABCD_subj_ls))
    ABCD_subj_ls = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt');
end

%% read text files
[ABCD_bhvrs, nbhvr] = CBIG_text2cell(ABCD_bhvr_ls);
HCP_bhvrs = CBIG_text2cell(HCP_bhvr_ls);
assert(length(HCP_bhvrs) == nbhvr, 'Number of behaviors in HCP_bhvr_ls and ABCD_bhvr_ls differ!')

ABCD_colloqs = CBIG_text2cell(ABCD_colloq_ls);
HCP_colloqs = CBIG_text2cell(HCP_colloq_ls);

d_abcd = readtable(ABCD_csv);
[ABCD_subj, nABCD] = CBIG_text2cell(ABCD_subj_ls);
[~, ~, idx] = intersect(ABCD_subj, d_abcd.subjectkey, 'stable');

alpha_FDR = 0.05;

metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
[~, midx] = intersect(metrics, test_metric, 'stable');
load(fullfile(train_dir, 'Pset.mat'))  % load multi-level block permutations
avg_stats = zeros(nbhvr, 1);
avg_null_stats = zeros(nbhvr, Nperm);
for b = 1:nbhvr
    ABCD_bhvr_nm = ABCD_bhvrs{b};   HCP_bhvr_nm = HCP_bhvrs{b};
    ABCD_y = d_abcd.(ABCD_bhvr_nm)(idx);

    curr_avg_stats = [];
    curr_avg_null_stats = [];
    for seed = 1:max_HCP_seed
        curr_train_dir = fullfile(train_dir, ['randseed_' num2str(seed)], HCP_bhvr_nm);
        curr_test_dir = fullfile(test_dir, ['randseed_' num2str(seed)], HCP_bhvr_nm);
        if(~exist(curr_train_dir, 'dir')); continue; end

        opt = load(fullfile(curr_train_dir, ['final_result_' HCP_bhvr_nm '.mat']));
        Nfolds = length(opt.optimal_acc);
        acc_out = fullfile(curr_test_dir, 'perm.mat');
        pred = load(fullfile(curr_test_dir, 'prediction_AAWA.mat'));
        curr_pred = [];
        for f = 1:Nfolds
            curr_pred = [curr_pred; pred.pred_stats{f}(midx)];
        end

        if(exist(acc_out, 'file'))
            load(acc_out)
        else
            for f = 1:Nfolds
                y_reg = load(fullfile(curr_train_dir, 'y', ['fold_' num2str(f)], ['y_regress_' HCP_bhvr_nm '.mat']));
                load(fullfile(curr_train_dir, 'FSM_innerloop', ['fold_' num2str(f)], 'FSM_corr.mat'))
                K_train = FSM;
                load(fullfile(testFSM_dir, ['randseed_' num2str(seed)], [HCP_bhvr_nm '_' ABCD_bhvr_nm], ...
                    ['fold_' num2str(f)], 'FSM_ABCD_full_vs_HCP.mat'))
                K_test = FSM;
                N_train = size(K_train, 1);
                N_test = size(K_test, 1);
                opt_lambda = opt.optimal_lambda(f);

                load(fullfile(curr_train_dir, ['no_relative_10_fold_sub_list_' HCP_bhvr_nm '.mat']))
                train_ind = sub_fold(f).fold_index==0;

                % compute the part of parameters that are not dependent on y
                % so that they are be shared across all permutations
                K_r = K_train + opt_lambda*eye(N_train);
                X = ones(N_train,1);
                inv_K_r = inv(K_r);
                beta_stable = (X' * (inv_K_r * X)) \ X' * inv_K_r;

                % for each permutation, calculate prediction accuracies
                for p = 2:Nperm+1
                    y_perm = y_reg.y_resid(Pset(:,p));
                    y_train = y_perm(train_ind);

                    beta = beta_stable * y_train;
                    alpha = inv_K_r * (y_train - X * beta);
            
                    y_p = K_test * alpha + ones(N_test,1) .* beta;
                    for k = 1:length(metrics)
                        stats_perm.(metrics{k})(f,p) = ...
                            CBIG_compute_prediction_acc_and_loss(y_p, ABCD_y, metrics{k}, y_train);
                    end
                end
            end
            save(acc_out, 'stats_perm')
        end
        curr_avg_stats = cat(1, curr_avg_stats, mean(curr_pred, 1));
        curr_avg_null_stats = cat(1, curr_avg_null_stats, mean(stats_perm.(test_metric)(:, 2:end), 1));
    end
    avg_stats(b) = mean(curr_avg_stats, 1);
    avg_null_stats(b,:) = mean(curr_avg_null_stats, 1);
    p_perm(b) = length(find(avg_null_stats(b,:) > avg_stats(b))) ./ Nperm;

    clear stats_perm
end

%% multiple comparisons correction
H_perm_FDR = FDR(p_perm, alpha_FDR);
sig_perm_idx = sort(H_perm_FDR);
sig_perm_HCPbhvr = HCP_bhvrs(sig_perm_idx);
sig_perm_HCPcolloq = HCP_colloqs(sig_perm_idx);
sig_perm_ABCDbhvr = ABCD_bhvrs(sig_perm_idx);
sig_perm_ABCDcolloq = ABCD_colloqs(sig_perm_idx);

outdir = fileparts(outmat);
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(outmat, 'alpha_FDR', 'p_perm', 'H_perm_FDR', 'sig_perm_idx', 'sig_perm_HCPbhvr', ...
    'sig_perm_HCPcolloq', 'sig_perm_ABCDbhvr', 'sig_perm_ABCDcolloq', 'avg_stats', 'avg_null_stats')

end
