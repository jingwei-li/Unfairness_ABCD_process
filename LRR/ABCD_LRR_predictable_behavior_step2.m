function ABCD_LRR_predictable_behavior_step2(model_dir, cov_stem, FC_file, folds_dir, folds_fstem, curr_fold, Nperm, bhvr_nm, colloq_nm)

% ABCD_LRR_predictable_behavior_step2(model_dir, FC_file, folds_dir, folds_fstem, Nperm, bhvr_nm, colloq_nm)
%
% Repeat linear ridge regression with permuted behavioral scores. This function needs 
% to be called for each behavioral measure. The permutations need to be pre-generated 
% using `ABCD_LRR_predictable_behavior_step1.m`.
%
% Inputs:
%   - model_dir
%     Directory of LRR results.
%   - FC_file
%     The RSFC mat file (full path) of all subjects involved in the cross-validation.
%   - folds_dir
%     Directory containing all mat files of fold splits.
%   - folds_fstem
%     Filename prefix of the fold splits mat files.
%   - Nperm
%     Number of permutations.
%   - bhvr_nm 
%     Current behavioral measure.
%
% Author: Jingwei Li

%% load data for LRR; calculate permuted KRR accuracies
load(FC_file)
features = reshape3d_to_2d(corr_mat);

%% Read covariates which need to be regressed out from RSFC
load(fullfile(model_dir, bhvr_nm, ['confounds_X' cov_stem '.mat']));

%% load the permutation indices
load(fullfile(model_dir, 'perm', 'Pset.mat'))

metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};

load(fullfile(folds_dir, ['sub_fold' folds_fstem '_' bhvr_nm '.mat']))
Nfolds = length(sub_fold);
    
for i = 1:length(metrics)
    %stats_perm.(metrics{i}) = zeros(Nfolds, Nperm);
    stats_perm.(metrics{i}) = zeros(1, Nperm);
end

%acc_out = fullfile(model_dir, 'perm', [bhvr_nm '.mat']);
%if(~exist(acc_out, 'file'))
    %for f = 1:Nfolds
    for f = curr_fold
        fprintf('fold: %d\n', f)
        perm_out = fullfile(model_dir, bhvr_nm, 'perm', ['fold_' num2str(f) '.mat']);
        if(exist(perm_out, 'file')); continue; end
        
        mkdir(fullfile(model_dir, bhvr_nm, 'perm'))
        test_ind = sub_fold(f).fold_index==1;
        train_ind = ~test_ind;

        y_reg = load(fullfile(model_dir, bhvr_nm, 'y', ['fold_' num2str(f)], ...
            ['y_regress_' bhvr_nm '.mat']));

        %% split training vs test data for RSFC
        feat_train = features(:, train_ind);
        feat_test = features(:, test_ind);

        %% do confound regression from features, if necessary
        if(~isempty(cov_X) && ~strcmpi(cov_X, 'none'))
            [feat_train, beta] = CBIG_regress_X_from_y_train(feat_train', ...
                cov_X(train_ind, :));
            beta_pre = load(fullfile(model_dir, bhvr_nm, ...
                'params', ['fold_' num2str(f)], 'feature_regress_beta.mat'));
            if(max(abs(beta(:) - beta_pre.beta(:))) > 1e-8)
                error('[Regression from RSFC]: beta differred from original elastic net results')
            end

            feat_test = CBIG_regress_X_from_y_test(feat_test', ...
                cov_X(test_ind, :), beta);
            feat_train = feat_train';  feat_test = feat_test';
        end

        opt_params = load(fullfile(model_dir, bhvr_nm, 'params', ['fold_' num2str(f)], ...
            ['selected_parameters_' bhvr_nm '.mat']));

        % for each permutation, calculate prediction accuracies
        for p = 2:Nperm+1
            rng('default')
            rng(1)

            % permute y
            y_perm = y_reg.y_resid(Pset(:,p));
            y_train = y_perm(train_ind);
            y_test = y_perm(test_ind);

            % select features
            if opt_params.curr_threshold ~= 1
                [feat_train, feat_test] = CBIG_FC_FeatSel( feat_train, feat_test, y_train, ...
                    opt_params.curr_threshold );
            end

            Mdl = fitrlinear(feat_train, y_train, 'ObservationsIn', 'columns', 'Lambda', ...
                opt_params.curr_lambda, 'Learner', 'leastsquares', 'Regularization','ridge');
            y_p = predict(Mdl, feat_test, 'ObservationsIn', 'columns');

            for k = 1:length(metrics)
                stats_perm.(metrics{k})(1,p) = ...
                    CBIG_compute_prediction_acc_and_loss(y_p, y_test, metrics{k}, y_train);
            end
        end
        save(perm_out, 'stats_perm')
    end
    %save(acc_out, 'stats_perm')
    %system(sprintf('ls -l %s', acc_out))
%end

end


function out = reshape3d_to_2d(features)
    % reshapes a #ROI x #ROI x #subjects matrix into
    % #ROI x #subjects by extracting the lower triangle
    % of the correlation
    temp = ones(size(features,1), size(features,2));
    tril_ind = tril(temp, -1);
    features_reshaped = reshape(features, size(features,1)*size(features,2), size(features, 3));
    out = features_reshaped(tril_ind==1, :);
end
