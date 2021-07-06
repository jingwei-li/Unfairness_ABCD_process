function ABCD_LRR_predictable_behavior_step3(model_dir, test_metric, outmat, bhvr_ls, colloq_ls, Nperm)

% ABCD_LRR_predictable_behavior_step3(model_dir, test_metric, outmat, bhvr_ls, colloq_ls)
%
% Test the true accuracy against the accuracy distributions after permutation. 
%
% Inputs:
%   - model_dir
%     Directory of LRR results.
%   - outmat
%     Output mat file containing the behaviors predicted more than chance.
%   - bhvr_ls  (optional)
%     List of behavioral measures.
%     Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt
%   - colloq_ls (optional)
%     List of colloquial names.
%     Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt
%
% Author: Jingwei Li

%% setup default parameters
ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
if(~exist('Nperm', 'var') || isempty(Nperm))
    Nperm = 1000;
end
alpha_FDR = 0.05;

%% read lists
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
[colloq_nm, ncolloq] = CBIG_text2cell(colloq_ls);
if(nbhvr ~= ncolloq)
    error('Number of behavioral names and number of colloquial names are not equal.')
end

%% compare the empirical accuracy and accuracy after permutation
p_perm = zeros(nbhvr, 1);
for b = 1:nbhvr
    opt = load(fullfile(model_dir, bhvr_nm{b}, 'results', 'optimal_acc', ...
        [bhvr_nm{b} '_final_acc.mat']));
    Nfolds = length(opt.optimal_statistics);
    orig_stats = zeros(Nfolds, 1);
    for f = 1:Nfolds
        orig_stats(f) = opt.optimal_statistics{f}.(test_metric);
    end

    acc_out = fullfile(model_dir, 'perm', [bhvr_nm{b} '.mat']);
    if(exist(acc_out, 'file'))
        load(acc_out)
    else
        metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
        for k = 1:length(metrics)
            stats_perm.(metrics{k}) = zeros(Nfolds, Nperm);
        end
        for f = 1:Nfolds
            curr_perm = load(fullfile(model_dir, bhvr_nm, 'perm', ['fold_' num2str(f) '.mat']));
            for k = 1:length(metrics)
                stats_perm.(metrics{k})(f,:) = curr_perm.stats_perm.(metrics{k})(1,:);
            end
        end
        save(acc_out, 'stats_perm')
    end

    avg_stats = mean(orig_stats, 1);
    avg_null_stats = mean(stats_perm.(test_metric)(:, 2:end), 1);
    p_perm(b) = length(find(avg_null_stats > avg_stats)) ./ Nperm;
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