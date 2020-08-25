function ABCD_KRR_corr_allAA_vs_randWA(model_dir, bhvr_ls, subj_ls, split_dir, split_fstem, ...
    Nsplits, outmat)

%% default arguments
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, '/behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

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

%% compute Pearson's correlation
corr_AA = nan(nbhvr, Nsplits); corr_WA = corr_AA; corr_AAWA = corr_AA;
AA_pred = cell(nbhvr, Nsplits); WA_pred = AA_pred;
AA_test = AA_pred; WA_test = AA_pred;
for b = 1:nbhvr
    fprintf('#%d behavior: %s ...\n', b, bhvr_nm{b})
    fold_file = fullfile(split_dir, 'allAA_randWA', [bhvr_nm{b} split_fstem '.mat']);
    load(fold_file)
    if(length(sub_fold) ~= Nsplits)
        error('Nsplits does not equal to length of sub_fold.')
    end
    
    opt_file = fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']);
    opt = load(opt_file);

    for f = 1:Nsplits
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

end