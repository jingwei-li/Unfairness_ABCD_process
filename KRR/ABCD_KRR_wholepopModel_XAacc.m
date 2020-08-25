function ABCD_KRR_wholepopModel_XAacc(XA, model_dir, split_dir, split_fstem, outstem, csvname, bhvr_ls)

% Syntax: ABCD_KRR_wholepopModel_XAacc(XA, model_dir, bhvr_ls, outmat)
%
% 

metrics = {'corr', 'COD', 'predictive_COD', 'MAE', 'MAE_norm', 'MSE', 'MSE_norm'};
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('cvsname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phtnotypes_pass_rs.txt');
end
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

d = readtable(csvname);
races = d.race;
switch XA
case 'WA'
    raceXA = 1;
case 'AA'
    raceXA = 2;
otherwise
    error('Not supporting race %s', XA)
end

for b = 1:nbhvr
    load(fullfile(model_dir, ['sub_fold' split_fstem bhvr_nm{b} '.mat']))
    opt = load(fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']));
    Nsplits = length(sub_fold);
    % initialize optimal stats structure
    for m = 1:length(metrics)
        optimal_stats.(metrics{m}) = zeros(Nsplits, 1);
    end

    for f = 1:Nsplits
        % load regressed y
        krry = load(fullfile(model_dir, 'y', ['fold_' num2str(f)], ...
            ['y_regress_' bhvr_n
        % load test CV results
        testcv = load(fullfile(model_dir, 'test_cv', ['fold_' num2str(f)], ...
            ['acc_' bhvr_nm{b} '.mat']));

        % race of test subjects of current fold
        [~, ~, idx] = intersect(sub_fold(f).subject_list, d.subjectkey, 'stable');
        curr_races = races(idx);
        XAidx = curr_races == raceXA;

        % collect true & predicted scores of test XA subjects
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
        curr_y_pred = testcv.y_p{opt_kernel_idx, opt_lambda_idx, opt_thres_idx}{1}(XAidx);
        curr_y_true = testcv.y_t{opt_kernel_idx, opt_lambda_idx, opt_thres_idx}{1}(XAidx);

        % collect training scores
        curr_y_train = krry.y_resid(sub_fold(f).fold_index==0);

        % compute stats
        for m = 1:length(metrics)
            optimal_stats.(metrics{m})(f) = CBIG_compute_prediction_acc_and_loss(curr_y_pred, curr_y_test, metrics{m}, curr_y_train);
        end
    end
    save(fullfile(model_dir, ['final_result' outstem '_' bhvr_nm{b} '.mat']))
    clear sub_fold opt
end
    
end