function [yp_grp2, yt_grp2, acc, pred_stats, y_grp2, y_grp2_resid, y_grp1_train_resid, y_grp3_resid] = ABCD_KRR_test_grp1Model_on_grp2_perfold( ...
	f, model_dir, behavior, split_grp1, idx_grp1_all, idx_grp2, cov_grp1, all_cov, all_y, corr_mat, opt, metrics, idx_grp3)

% [yp_grp2, yt_grp2, acc, pred_stats, y_grp2_resid] = ABCD_KRR_test_grp1Model_on_grp2_perfold( ...
%     f, model_dir, behavior, split_grp1, idx_grp1_all, idx_grp2, cov_grp1, all_cov, all_y, corr_mat, opt, metrics)
%
% Inputs:
%   - f: scalar. Fold index.
%   - model_dir: directory to the models trained on grp1 (full path). It contains a subfolder for each behavior.
%   - behavior: behavior name.
%   - split_grp1: a structure with lengths of total number of folds/splits. The splits of grp1 subjects. 
%                 It contains a field fold_index, where split_grp1(f).fold_index == 0 means the corresponding 
%                 subject was in the training set.
%   - idx_grp1_all: indices of grp1 subjects in all subjects used for the whole project.
%   - idx_grp2: indices of grp2 subjects in all subjects used for the whole project
%   - cov_grp1: #grp1 subjects x #confounds matrix. Confounding variables of grp1 subjects.
%   - all_cov: #all subjects x #confounds matrix. Confounding variables of all subjects used for the 
%              whole project.
%   - all_y: column vector with length = #all subjects. Raw behavioral scores of all subjects used 
%            for the whole project.
%   - corr_mat: (#ROI * (#ROI-1) / 2) x #all subjects matrix. Vectorized functional connectivity of
%               all subjects used for the whole project.
%   - opt: optimal results structure. Loaded from 
%          fullfile(model_dir, behavior, ['final_result_' behavior '.mat'])
%   - metrics: a cell containing the prediction statistics. 
%              E.g. {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'}

y_fold_grp1 = load(fullfile(model_dir, behavior, 'y', ['fold_' num2str(f)], ...
	['y_regress_' behavior '.mat']));
y_grp1_train = y_fold_grp1.y_orig(split_grp1(f).fold_index == 0);
y_grp2 = all_y(idx_grp2);
if(exist('idx_grp3', 'var') && ~isempty(idx_grp3))
	y_grp3 = all_y(idx_grp3);
end

%% regression, if necessary
y_grp3_resid = [];
if(strcmpi(cov_grp1, 'none'))
	y_grp2_resid = y_grp2;
	y_grp1_train_resid = y_grp1_train;
	if(exist('idx_grp3', 'var') && ~isempty(idx_grp3))
		y_grp3_resid = y_grp3;
	end
else
	if(isempty(cov_grp1))
		cov_grp1_train = [];
		cov_grp2 = [];
		if(exist('idx_grp3', 'var') && ~isempty(idx_grp3))
			cov_grp3 = [];
		end
	else
		cov_grp1_train = cov_grp1(split_grp1(f).fold_index == 0, :);
		cov_grp2 = all_cov(idx_grp2, :);
		if(exist('idx_grp3', 'var') && ~isempty(idx_grp3))
			cov_grp3 = all_cov(idx_grp3, :);
		end
	end

	[y_grp1_train_resid, beta] = CBIG_regress_X_from_y_train(y_grp1_train, cov_grp1_train);
	if(max(abs(y_grp1_train_resid - y_fold_grp1.y_resid(split_grp1(f).fold_index == 0))) > 1e-6);
		error('Regression in training sample was not replicated. Max diff: %f', max(abs(y_grp1_train_resid - y_fold_grp1.y_resid(split_grp1(f).fold_index == 0))) )
	end
	y_grp2_resid = CBIG_regress_X_from_y_test(y_grp2, cov_grp2, beta);
	if(exist('idx_grp3', 'var') && ~isempty(idx_grp3))
		y_grp3_resid = CBIG_regress_X_from_y_test(y_grp3, cov_grp3, beta);
	end
end

%% compute kernel between training subjects in grp1 and testing subjects in grp2.
idx_grp1_train = idx_grp1_all(split_grp1(f).fold_index == 0);
FC_grp1_train = corr_mat(:, idx_grp1_train);
FC_grp2 = corr_mat(:, idx_grp2);
FSM = CBIG_crossvalid_kernel_with_scale(FC_grp1_train, FC_grp2, ...
	[], [], opt.optimal_kernel(f).type, opt.optimal_kernel(f).scale);
FSM = FSM( (length(idx_grp1_train)+1):end, 1:length(idx_grp1_train) );

%% load kernel among training subjects in grp1
FSM_train_nm = ['FSM_' opt.optimal_kernel(f).type];
if(~strcmpi(opt.optimal_kernel(f).type, 'corr'))
	FSM_train_nm = [FSM_train_nm '_' num2str(opt.optimal_kernel(f).scale)];
end
FSM_train = load(fullfile(model_dir, behavior, 'FSM', [FSM_train_nm '.mat']));
FSM_train = FSM_train.FSM(split_grp1(f).fold_index == 0, split_grp1(f).fold_index == 0);

%% use optimal parameter to test grp2 subjects
[yp_grp2, yt_grp2, acc, pred_stats] = CBIG_KRR_test_cv( 0, FSM_train, FSM, ...
	y_grp1_train_resid, y_grp2_resid, y_grp2, 1, ...
	opt.optimal_lambda(f), opt.optimal_threshold(f), metrics );
	
end