function ABCD_KRR_errors_AAvsWA(model_dir, bhvr_ls, colloq_ls, subj_ls, split_dir, split_fstem, Nsplits, ...
    metric, outmat, bin_flag)

% ABCD_KRR_errors_AAvsWA(model_dir, bhvr_ls, colloq_ls, subj_ls, split_dir, split_fstem, Nsplits, ...
%     metric, outmat, bin_flag)
%
% 

%% default arguments
ls_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';

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

%% compute predictive COD
pCOD_AA = nan(nbhvr, Nsplits); pCOD_WA = pCOD_AA; pCOD_AAWA = pCOD_AA;
ss_res_AA = pCOD_AA; ss_res_WA = pCOD_WA; ss_res_AAWA = pCOD_AA;
ss_total = pCOD_AA;

AA_pred = cell(nbhvr, Nsplits); WA_pred = AA_pred; AA_test = AA_pred;
WA_test = AA_pred; AA_train = AA_pred; WA_train = AA_pred;
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
        
        %% collect true scores of training AA or WA subjects
        trainAA = setdiff(all_selAA{b}, sub_fold(f).selAA);
        trainWA = setdiff(all_selWA{b}, sub_fold(f).selWA);
        
        [~, idx] = intersect(subjects, trainAA, 'stable');
        AA_train{b,f} = krry.y_resid(idx);
        [~, idx] = intersect(subjects, trainWA, 'stable');
        WA_train{b,f} = krry.y_resid(idx);
        
        %% compute predictive COD
        [err_AA(b,f), err_WA(b,f), err_AAWA(b,f)] = ...
            ABCD_errors_2groups(AA_pred{b,f}, WA_pred{b,f}, AA_test{b,f}, ...
            WA_test{b,f}, metric, AA_train{b,f}, WA_train{b,f});
    end
end

switch metric
case 'MSE'
    MSE_AA = err_AA;  MSE_WA = err_WA;  MSE_AAWA = err_AAWA;
case 'MSE_norm'
    MSE_norm_AA = err_AA;  MSE_norm_WA = err_WA;  MSE_norm_AAWA = err_AAWA;
case 'MAE'
    MAE_AA = err_AA;  MAE_WA = err_WA;  MAE_AAWA = err_AAWA;
case 'MAE_norm'
    MAE_norm_AA = err_AA;  MAE_norm_WA = err_WA;  MAE_norm_AAWA = err_AAWA;
end
save(outmat, [metric '_AA'], [metric '_WA'], [metric '_AAWA'], ...
    'ss_total', 'AA_pred', 'WA_pred', 'AA_test', 'WA_test', 'AA_train', 'WA_train');
    
end