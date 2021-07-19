function ABCD_KRR_test_HCP_predictable_behavior( outmat, ABCD_split_dir, train_dir, test_dir, testFSM_dir, Nperm, test_metric, ...
    ABCD_bhvr_ls, ABCD_colloq_ls, HCP_bhvr_ls, HCP_colloq_ls, HCP_KRR_dir )

% ABCD_KRR_test_HCP_predictable_behavior( outmat, train_dir, test_dir, testFSM_dir, Nperm, test_metric, ...
%     ABCD_bhvr_ls, ABCD_colloq_ls, HCP_bhvr_ls, HCP_colloq_ls, HCP_KRR_dir )
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

%% read text files
[ABCD_bhvrs, nbhvr] = CBIG_text2cell(ABCD_bhvr_ls);
HCP_bhvrs = CBIG_text2cell(HCP_bhvr_ls);
assert(length(HCP_bhvrs) == nbhvr, 'Number of behaviors in HCP_bhvr_ls and ABCD_bhvr_ls differ!')
ABCD_colloqs = CBIG_text2cell(ABCD_colloq_ls);
HCP_colloqs = CBIG_text2cell(HCP_colloq_ls);

alpha_FDR = 0.05;
metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
[~, midx] = intersect(metrics, test_metric);
load(fullfile(train_dir, 'perm', 'Pset.mat'))  % load multi-level block permutations
avg_stats = zeros(nbhvr, 1);
avg_null_stats = zeros(nbhvr, Nperm);
for b = 1:nbhvr
    ABCD_bhvr_nm = ABCD_bhvrs{b};   HCP_bhvr_nm = HCP_bhvrs{b};
    HCP_y = load(fullfile(HCP_KRR_dir, ['y_' HCP_bhvr_nm '.mat']));

    opt = load(fullfile(train_dir, ABCD_bhvr_nm, ['final_result_' ABCD_bhvr_nm '.mat']));
    Nfolds = length(opt.optimal_acc);
    acc_out = fullfile(test_dir, 'perm.mat');
    if(exist(acc_out, 'file'))
        load(acc_out)
    else
        for f = 1:Nfolds
            y_reg = load(fullfile(train_dir, ABCD_bhvr_nm, 'y', ['fold_' num2str(f)], ['y_regress_' ABCD_bhvr_nm '.mat']));
            load(fullfile(train_dir, ABCD_bhvr_nm, 'FSM_innerloop', ['fold_' num2str(f)], 'FSM_corr.mat'))
            K_train = FSM;
            load(fullfile(testFSM_dir, [ABCD_bhvr_nm '_' HCP_bhvr_nm], ['fold_' num2str(f)], 'FSM_HCP_full_vs_ABCD.mat'))

            K_test = FSM;
            N_train = size(K_train, 1);
            N_test = size(K_test, 1);
            opt_lambda = opt.optimal_lambda(f);

            load(fullfile(ABCD_split_dir, ['sub_fold_pass_rs_pass_pheno_' ABCD_bhvr_nm '.mat']))
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
                        CBIG_compute_prediction_acc_and_loss(y_p, HCP_y.y, metrics{k}, y_train);
                end
            end
        end
        save(acc_out, 'stats_perm')
    end

    pred = load(fullfile(test_dir, ABCD_bhvr_nm, 'prediction.mat'));
    curr_pred = [];
    for f = 1:Nfolds
        curr_pred = [curr_pred; pred.pred_stats{f}(midx)];
    end
    avg_stats(b) = mean(curr_pred, 1);
    avg_null_stats(b, :) = mean(stats_perm.(test_metric)(:, 2:end), 1);
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