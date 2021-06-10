function ABCD_KRR_predictable_behavior(model_dir, folds_dir, folds_fstem, Nperm, test_metric, outmat, bhvr_ls, colloq_ls, subj_ls, csvname)

% ABCD_KRR_predictable_behavior(model_dir, cov_stem, folds_dir, folds_fstem, outmat, bhvr_ls, subj_ls, csvname)
%
% Permute behavioral scores across subjects within each site, to test the predictability
% of the original kernal regression model.
%
% Inputs:
%   - model_dir
%     Directory of KRR results.
%   - folds_dir
%     Directory containing all mat files of fold splits.
%   - folds_fstem
%     Filename prefix of the fold splits mat files.
%   - Nperm
%     Number of permutations.
%   - test_metric
%     The metric used to perform statistical testing.
%   - outmat
%     Output mat file containing the behaviors predicted more than chance.
%   - bhvr_ls  (optional)
%     List of behavioral measures.
%     Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt
%   - colloq_ls (optional)
%     List of colloquial names.
%     Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt
%   - subj_ls  (optional)
%     Subject list.
%     Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt
%   - csvname  (optional)
%     Name of the csv file containing site information.
%     Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/phenotypes_pass_rs.txt
%

%% setup default parameters
ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end
alpha_FDR = 0.05;

%% read shared files
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
[colloq_nm, ncolloq] = CBIG_text2cell(colloq_ls);
if(nbhvr ~= ncolloq)
    error('Number of behavioral names and number of colloquial names are not equal.')
end
[subjects, nsub] = CBIG_text2cell(subj_ls);
d = readtable(csvname);
subj_hdr = 'subjectkey';
site_hdr = 'site';
fam_hdr = 'family_id';
grp_hdr = 'rel_grp_id';

[~, ~, idx] = intersect(subjects, d.(subj_hdr), 'stable');
site = d.(site_hdr)(idx);
fam = d.(fam_hdr)(idx);
grp = d.(grp_hdr)(idx);

mkdir(fullfile(model_dir, 'perm'))

%% generate exchangeable blocks and permutations
idx_perm_out = fullfile(model_dir, 'perm', 'Pset.mat');
if(~exist(idx_perm_out, 'file'))
    [Pset, B] = ABCD_multilevel_block_perm(site, fam, grp, Nperm+1);
    N_uq_perm = size(unique(Pset', 'rows'), 1);
    if (N_uq_perm-1 < Nperm)
        warning('Number of required permutations is higher than maximal possible permutations.')
    end
    save(idx_perm_out, 'Pset', 'B')
else
    load(idx_perm_out)
end

%% load data for KRR; calculate permuted KRR accuracies
metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
p_perm = zeros(nbhvr, 1);
[flag, msg] = system(['ls ' fullfile(model_dir, 'final_result*.mat')])
for b = 1:nbhvr
    load(fullfile(folds_dir, ['sub_fold' folds_fstem '_' bhvr_nm{b} '.mat']))
    Nfolds = length(sub_fold);
    if(~flag)
        opt = load(fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']));
    else
        opt = load(fullfile(model_dir, bhvr_nm{b}, ['final_result_' bhvr_nm{b}, '.mat']));
    end

    for i = 1:length(metrics)
        stats_perm.(metrics{i}) = zeros(Nfolds, Nperm);
    end

    acc_out = fullfile(model_dir, 'perm', [bhvr_nm{b} '.mat']);
    if(exist(acc_out, 'file'))
        load(acc_out)
    else
        for f = 1:Nfolds
            test_ind = sub_fold(f).fold_index==1;
            train_ind = ~test_ind;
            N_train = length(find(train_ind));
            N_test = length(find(test_ind));

            if(~flag)
                y_reg = load(fullfile(model_dir, 'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']));
            else
                y_reg = load(fullfile(model_dir, bhvr_nm{b}, 'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']));
            end

            if(~flag)
                load(fullfile(model_dir, 'FSM', 'FSM_corr.mat'))
                K_train = FSM(train_ind, train_ind);
                K_test = FSM(test_ind, train_ind);
                clear FSM
            else
                load(fullfile(model_dir, bhvr_nm{b}, 'FSM_innerloop', ['fold_' num2str(f)], 'FSM_corr.mat'))
                K_train = FSM;
                clear FSM
                load(fullfile(model_dir, bhvr_nm{b}, 'FSM_test', ['fold_' num2str(f)], 'FSM_corr.mat'))
                K_test = FSM(test_ind, train_ind);
                clear FSM
            end

            opt_lambda = opt.optimal_lambda(f);

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
                y_test = y_perm(test_ind);

                beta = beta_stable * y_train;
                alpha = inv_K_r * (y_train - X * beta);
            
                y_p = K_test * alpha + ones(N_test,1) .* beta;
                for k = 1:length(metrics)
                    stats_perm.(metrics{k})(f,p) = ...
                        CBIG_compute_prediction_acc_and_loss(y_p, y_test, metrics{k}, y_train);
                end
            end

        end
        save(acc_out, 'stats_perm')

    end

    avg_stats = mean(opt.optimal_stats.(test_metric), 1);
    avg_null_stats = mean(stats_perm.(test_metric)(:, 2:end), 1);
    p_perm(b) = length(find(avg_null_stats > avg_stats)) ./ Nperm;

    clear stats_perm
end

%% multiple comparisons correction
H_perm_FDR = FDR(p_perm, alpha_FDR);
sig_perm_idx = sort(H_perm_FDR);
sig_perm_bhvr = bhvr_nm(sig_perm_idx);
sig_perm_colloq = colloq_nm(sig_perm_idx);

outdir = fileparts(outmat);
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(outmat, 'alpha_FDR', 'p_perm', 'H_perm_FDR', 'sig_perm_idx', 'sig_perm_bhvr', ...
    'sig_perm_colloq')
    
end




