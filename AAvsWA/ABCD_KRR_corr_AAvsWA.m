function ABCD_KRR_corr_AAvsWA(model_dir, bhvr_ls, colloq_ls, subj_ls, split_dir, split_fstem, Nsplits, ...
    predictable_stats, r_th_ls, outmat)

% ABCD_KRR_corr_AAvsWA(model_dir, bhvr_ls, colloq_ls, subj_ls, split_dir, split_fstem, Nsplits, ...
%     predictable_stats, outmat)
% 
% Calculate Pearson's correlation between predicted behavioral scores and true behavioral scores of
% matched AA and WA for the kernel ridge regression models trained on whole population.
%
% Input:
%   - model_dir
%     The directory storing kernel regression results (full path).
%
%   - bhvr_ls (optional)
%     Behavior list (full path, text file).
%     Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
% 
%   - colloq_ls (optional)
%     List of behaviors' colloquial name (full path).
%     Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt'
%
%   - subj_ls (optional)
%     Subject list (full path). Default: 
%  '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt'
%
%   - split_dir
%     The directory storing the data split (full path). It contains a file
%     called ['sel_AAWA' split_fstem '.mat']. For each behaviors, there
%     should be a file called ['sub_fold' split_fstem '_' <behavior name>
%     '.mat'].
%
%   - split_fstem
%     The string that was attached to the filenames of the data split
%     files.
%
%   - Nsplits (optional)
%     Number of splits. Default: 120 (10 choose 3).
%  
%   - predictable_stats
%     Filename of the permutation test result of predictability (full
%     path).
%
%   - r_th_ls
%     Filename of the list of behaviors whose correlation accuracy across all participants is above a threshold.
%
%   - outmat
%     Output filename (full path).
%
% Author: Jingwei Li

%% default arguments
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, '/behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
colloq_nm = CBIG_text2cell(colloq_ls);

if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_ls);

if(~exist('Nsplits', 'var') || isempty(Nsplits))
    Nsplits = nchoosek(10, 3);
end

%% set default hyperparamters if not passed in
if(~exist('ker_param', 'var') || strcmpi(ker_param, 'none'))
    ker_param.type = 'corr';
    ker_param.scale = NaN;
end
ker_param = struct2cell(ker_param);

if(~exist('lambda_set', 'var') || strcmpi(lambda_set, 'none'))
    lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
        5 10 15 20];
end

if(~exist('bin_flag', 'var') || isempty(bin_flag))
    bin_flag = 0;
end
if(bin_flag==1)
    if(~exist('threshold_set', 'var') || strcmpi(threshold_set, 'none') || isempty(threshold_set))
        threshold_set = [-1:0.1:1];
    end
else
    threshold_set = NaN;
end

%% load all seletect AA, WA
match = load(fullfile(split_dir, ['sel_AAWA' split_fstem '.mat']));
all_selAA = cell(nbhvr,1); all_selWA = all_selAA;
for b = 1:nbhvr
    all_selAA{b} = cat(2, match.selAA{b,:});
    all_selWA{b} = cat(1, match.selWA{b,:});
end

%% compute Pearson's correlation
corr_AA = nan(nbhvr, Nsplits); corr_WA = corr_AA; corr_AAWA = corr_AA;
AA_pred = cell(nbhvr, Nsplits); WA_pred = AA_pred;
AA_test = AA_pred; WA_test = AA_pred;
for b = 1:nbhvr
    fprintf('#%d behavior: %s ...\n', b, bhvr_nm{b})
    fold_file = fullfile(split_dir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']);
    load(fold_file)
    if(length(sub_fold) ~= Nsplits)
        error('Nsplits does not equal to length of sub_fold.')
    end

    opt_file = fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']);
    opt = load(opt_file);

    for f = 1:length(sub_fold)
        krry = load(fullfile(model_dir, 'y', ['fold_' num2str(f)], ...
            ['y_regress_' bhvr_nm{b} '.mat']));
        testcv = load(fullfile(model_dir, 'test_cv', ['fold_' num2str(f)], ...
            ['acc_' bhvr_nm{b} '.mat']));

        %% collect true & predicted scores of test AA or WA subjects
        AAidx = zeros(length(sub_fold(f).subject_list), 1);
        WAidx = zeros(length(sub_fold(f).subject_list), 1);
        [~, idx] = intersect(sub_fold(f).subject_list, sub_fold(f).selAA, 'stable');
        AAidx(idx) = 1;
        [~, idx] = intersect(sub_fold(f).subject_list, sub_fold(f).selWA, 'stable');
        WAidx(idx) = 1;
        
        if(strcmp(opt.optimal_kernel(f).type, 'corr'))
            opt_kernel_idx = strcmp(ker_param(1,:,:), opt.optimal_kernel(f).type);
        else
            opt_kernel_idx = strcmp(ker_param(1,:,:), opt.optimal_kernel(f).type) ...
                & cell2mat(ker_param(2,:,:)) == opt.optimal_kernel(f).scale;
        end
        opt_lambda_idx = lambda_set == opt.optimal_lambda(f);
        if(bin_flag==1)
            opt_thres_idx = threshold_set == opt.optimal_threshold(f);
        else
            opt_thres_idx = 1;
        end
        AA_pred{b,f} = testcv.y_p{opt_kernel_idx, opt_lambda_idx, opt_thres_idx}{1}(logical(AAidx));
        WA_pred{b,f} = testcv.y_p{opt_kernel_idx, opt_lambda_idx, opt_thres_idx}{1}(logical(WAidx));
        
        AA_test{b,f} = testcv.y_t{opt_kernel_idx, opt_lambda_idx, opt_thres_idx}{1}(logical(AAidx));
        WA_test{b,f} = testcv.y_t{opt_kernel_idx, opt_lambda_idx, opt_thres_idx}{1}(logical(WAidx));
        
        %% compute Pearson's r
        corr_AA(b,f) = CBIG_corr(AA_pred{b,f}, AA_test{b,f});
        corr_WA(b,f) = CBIG_corr(WA_pred{b,f}, WA_test{b,f});
        corr_AAWA(b,f) = CBIG_corr([AA_pred{b,f}; WA_pred{b,f}], ...
            [AA_test{b,f}; WA_test{b,f}]);
    end
end

save(outmat, 'corr_AA', 'corr_WA', 'corr_AAWA', 'AA_pred', 'WA_pred', ...
    'AA_test', 'WA_test')

%% select behaviors with significant Pearson's r, 
% and positive in either AA or WA
% and whose correlation accuracy across all participants is above a threshold (e.g. 0.15)
stats = load(predictable_stats);
sig_idx = zeros(nbhvr,1);
sig_idx(stats.sig_perm_idx) = 1;
idx = sig_idx==1 & (mean(corr_AA,2)>0 | mean(corr_WA,2)>0);
predictable = bhvr_nm(idx);
r_th_bhvr = CBIG_text2cell(r_th_ls);
[~, rth_idx] = intersect(predictable, r_th_bhvr, 'stable');
predictable = predictable(rth_idx);
if(~exist(fullfile(model_dir, 'lists'), 'dir'))
    mkdir(fullfile(model_dir, 'lists'))
end
CBIG_cell2text(predictable, fullfile(model_dir, 'lists', 'corr_predictable.txt'))
predictable = colloq_nm(idx);
predictable = predictable(rth_idx);
CBIG_cell2text(predictable, fullfile(model_dir, 'lists', 'corr_predictable_colloquial.txt'))

[outdir, outname, outext] = fileparts(outmat);
corr_AA = corr_AA(idx,:);  corr_WA = corr_WA(idx,:); corr_AAWA = corr_AAWA(idx,:);
AA_pred = AA_pred(idx,:);  WA_pred = WA_pred(idx,:);
AA_test = AA_test(idx,:);  WA_test = WA_test(idx,:);

corr_AA = corr_AA(rth_idx,:);  corr_WA = corr_WA(rth_idx,:); corr_AAWA = corr_AAWA(rth_idx,:);
AA_pred = AA_pred(rth_idx,:);  WA_pred = WA_pred(rth_idx,:);
AA_test = AA_test(rth_idx,:);  WA_test = WA_test(rth_idx,:);
save(fullfile(outdir, [outname '_predictable' outext]), ...
    'corr_AA', 'corr_WA', 'corr_AAWA', 'AA_pred', 'WA_pred', 'AA_test', 'WA_test')

end

