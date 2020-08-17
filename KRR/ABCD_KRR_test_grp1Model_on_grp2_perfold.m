function [yp_grp2, yt_grp2, acc, pred_stats, y_grp2, y_grp2_resid] = ABCD_KRR_test_grp1Model_on_grp2_perfold( ...
	f, model_dir, behavior, split_grp1, idx_grp1_all, idx_grp2, cov_grp1, all_cov, all_y, corr_mat, opt, metrics)

% [yp_grp2, yt_grp2, acc, pred_stats, y_grp2_resid] = ABCD_KRR_test_grp1Model_on_grp2_perfold( ...
%     split_grp1, y_fold_grp1, idx_grp1_all, idx_grp2, cov_grp1, all_cov, corr_mat, opt)
%
% Long description

y_fold_grp1 = load(fullfile(model_dir, behavior, 'y', ['fold_' num2str(f)], ...
	['y_regress_' behavior '.mat']));
y_grp1_train = y_fold_grp1.y_orig(split_grp1(f).fold_index == 0);
y_grp2 = all_y(idx_grp2);

%% regression, if necessary
if(strcmpi(cov_grp1, 'none'))
	y_grp2_resid = y_grp2;
else
	if(isempty(cov_grp1))
		cov_grp1_train = [];
		cov_grp2 = [];
	else
		cov_grp1_train = cov_grp1(split_grp1(f).fold_index == 0, :);
		cov_grp2 = all_cov(idx_grp2, :);
	end

	[y_grp1_train_resid, beta] = CBIG_regress_X_from_y_train(y_grp1_train, cov_grp1_train);
	if(~isequal(y_grp1_train_resid, y_fold_grp1.y_resid(split_grp1(f).fold_index == 0)));
		error('Regression in training sample was not replicated.')
	end
	y_grp2_resid = CBIG_regress_X_from_y_test(y_grp2, cov_grp2, beta);
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