function ABCD_KRR_wholepopModel_XAacc(XA, model_dir, split_dir, split_fstem, outstem, csvname, bhvr_ls)

% Syntax: ABCD_KRR_wholepopModel_XAacc(XA, model_dir, bhvr_ls, outmat)
%
% 

metrics = {'corr', 'COD', 'predictive_COD', 'MAE', 'MAE_norm', 'MSE', 'MSE_norm'};
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('cvsname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

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

%% read csv file, configure race of XA
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
    fprintf('#%d behavior: %s\n', b, bhvr_nm{b})
    load(fullfile(split_dir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']))
    opt = load(fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']));
    Nsplits = length(sub_fold);
    % initialize optimal stats structure
    for m = 1:length(metrics)
        optimal_stats.(metrics{m}) = zeros(Nsplits, 1);
    end

    y_pred = cell(Nsplits, 1); y_true = y_pred; y_train = y_pred;
    for f = 1:Nsplits
        % load regressed y
        krry = load(fullfile(model_dir, 'y', ['fold_' num2str(f)], ...
            ['y_regress_' bhvr_nm{b} '.mat']));
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
            optimal_stats.(metrics{m})(f) = CBIG_compute_prediction_acc_and_loss(curr_y_pred, curr_y_true, metrics{m}, curr_y_train);
        end
        y_pred{f} = curr_y_pred; y_true{f} = curr_y_true; y_train{f} = curr_y_train;
    end
    save(fullfile(model_dir, ['final_result' outstem '_' bhvr_nm{b} '.mat']), ...
        'optimal_stats', 'y_pred', 'y_true', 'y_train')
    clear sub_fold opt y_pred y_true y_train
end
    
end