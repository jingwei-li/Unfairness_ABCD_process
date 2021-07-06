function ABCD_KRR_training_prediction(model_dir, has_subdir, split_dir, split_fstem, bhvr_ls, LITE)

% ABCD_KRR_training_prediction(model_dir, has_subdir, split_dir, split_fstem, bhvr_ls)
%
% Obtain predicted score in the training set with optimal hyperparameters.

%% KRR default parameters
with_bias = 1;
bin_flag= 0;
metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};

%% Default of optional arguments
if(ischar(has_subdir))
    has_subdir = str2num(has_subdir);
end
if(~exist('LITE', 'var') || isempty(LITE))
    LITE = true;
end

proj_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race');
ls_dir = fullfile(proj_dir, 'scripts', 'lists');
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

%% load FSM, if KRR results WERE NOT saved for each behavior separately.
if(~has_subdir && LITE)
    load(fullfile(model_dir, 'FSM', 'FSM_corr.mat'))
end

for b = 1:nbhvr
    fprintf('#%d behavior: %s\n', b, bhvr_nm{b})
    %% load fold splits
    split_fname = fullfile(split_dir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']);
    if(~exist(split_fname, 'file'))
        split_fname = fullfile(split_dir, [bhvr_nm{b} split_fstem '.mat']);
    end
    load(split_fname)
    nfolds = length(sub_fold);

    %% load KRR optimal hyperparameters
    if(has_subdir)
        opt = fullfile(model_dir, bhvr_nm{b}, ['final_result_' bhvr_nm{b} '.mat']);
    else
        opt = fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']);
    end
    opt = load(opt);

    %% load FSM, if KRR results WERE saved for each behavior separately.
    if(has_subdir && LITE)
        load(fullfile(model_dir, bhvr_nm{b}, 'FSM', 'FSM_corr.mat'))
    end

    for f = 1:nfolds
        fprintf('\tFold %d\n', f)
        if(~LITE)
            if(has_subdir)
                load(fullfile(model_dir, bhvr_nm{b}, 'FSM_innerloop', ['fold_' num2str(f)], 'FSM_corr.mat'))
            else
                load(fullfile(model_dir, 'FSM_innerloop', ['fold_' num2str(f)], 'FSM_corr.mat'))
            end
        end

        %% load regressed y, select regressed y of training subjects
        if(has_subdir)
            y = fullfile(model_dir, bhvr_nm{b}, 'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']);
        else
            y = fullfile(model_dir, 'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']);
        end
        y = load(y);
        y_resid = y.y_resid(sub_fold(f).fold_index==0);
        y_orig = y.y_orig(sub_fold(f).fold_index==0);

        %% select FSM among training subjects
        if(LITE)
            kernel_train = FSM(sub_fold(f).fold_index==0, sub_fold(f).fold_index==0);
        else
            kernel_train = FSM;
        end

        %% train KRR with optimal hyperparameters, and predict behaviors of training subjects
        [y_p, y_t, acc, pred_stats] = CBIG_KRR_test_cv_training_scores( bin_flag, kernel_train, ...
            y_resid, y_orig, with_bias, opt.optimal_lambda(f), opt.optimal_threshold(f), metrics );
        acc

        %% save prediction results of training subjects
        if(has_subdir)
            outdir = fullfile(model_dir, bhvr_nm{b}, 'test_cv', ['fold_' num2str(f)]);
        else
            outdir = fullfile(model_dir, 'test_cv', ['fold_' num2str(f)]);
        end
        save(fullfile(outdir, ['opt_training_set_' bhvr_nm{b}, '.mat']), 'y_p', 'y_t', 'acc', 'pred_stats')
    end
end
    
end