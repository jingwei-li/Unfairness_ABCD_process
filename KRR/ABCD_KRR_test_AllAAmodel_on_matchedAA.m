function ABCD_KRR_test_AllAAmodel_on_matchedAA(csvname, model_dir, bhvr_ls, cfds_ls, full_subj_ls, ...
	full_FC, split_dir, split_stem, AA_subdir)

% ABCD_KRR_test_AllAAmodel_on_matchedAA(csvname, model_dir, bhvr_ls, cfds_ls, full_subj_ls, ...
%     full_FC, split_dir, split_stem, AA_subdir)
%
% Long description

%% default input arguments
[csvname, subj_hdr, d, bhvr_nm, nbhvr, cfds_nm, ncfds, full_subj_ls, all_subj, idx_full_subj, corr_mat] = ...
	ABCD_KRR_test_grp1Model_on_grp2_parse_args(csvname, bhvr_ls, cfds_ls, full_subj_ls, full_FC);

%% load confounds of all subjects
fprintf('Loading confounds of all subjects...\n')
if(strcmpi(cfds_ls, 'none'))
	fprintf('No regressor to be regressed from behaviors.\n')
	covariates = 'NONE';
else
	cfds_types = repmat({'continuous'}, 1, ncfds);
	sex_idx = strcmpi(cfds_nm, 'sex');
	if(any(sex_idx))
		cfds_types{sex_idx} = 'categorical';
	end
	
	all_cov = CBIG_read_y_from_csv( {csvname}, subj_hdr, cfds_nm, cfds_types, full_subj_ls, 'NONE', ',' );
end

metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
for b = 1:nbhvr
	fprintf('#%d behavior: %s\n', b, bhvr_nm{b});
	split_matched = load(fullfile(split_dir, ['sub_fold' split_stem '_' bhvr_nm{b} '.mat']));
	split_allAA = load(fullfile(split_dir, AA_subdir, [bhvr_nm{b} split_stem '.mat']));
	allAA = CBIG_text2cell(fullfile(split_dir, AA_subdir, ['subj_' bhvr_nm{b} split_stem '.txt']));
	Nsplits = length(split_allAA.sub_fold);
	if(length(split_matched.sub_fold) ~= Nsplits)
		error('Lengths of split_matched and split_allAA are not equal.')
	end

	% load covariates of all AA
	load(fullfile(model_dir, bhvr_nm{b}, ['confounds_' strjoin(cfds_nm, '_') '.mat']));
	if(size(covariates, 2) ~= ncfds)
		error('Number of confounds in cfds_ls and model confounds .mat file not equal.')
	end

	% load optimal parameters
	opt = load(fullfile(model_dir, bhvr_nm{b}, ['final_result_' bhvr_nm{b} '.mat']));

	% collect raw behavioral scores of full set of subjects
	all_y = d.(bhvr_nm{b});
	all_y = all_y(idx_full_subj);

	[~, ~, idx_allAA_all] = intersect(allAA, all_subj, 'stable');
	yp_matchedAA = cell(Nsplits, 1);  yt_matchedAA = yp_matchedAA;
	y_matchedAA = yp_matchedAA;   y_matchedAA_resid = yp_matchedAA;
	optimal_acc = nan(Nsplits, 1);
    fprintf('Fold = ')
    for f = 1:Nsplits
        fprintf('%d..', f)
		[~, ~, idx_matchedAA] = intersect(split_matched.sub_fold(f).selAA, all_subj, 'stable');

		[yp, yt, optimal_acc(f), pred_stats, y_matchedAA{f}, y_matchedAA_resid{f}, y_train_resid{f}] = ...
			ABCD_KRR_test_grp1Model_on_grp2_perfold( f, model_dir, bhvr_nm{b}, split_allAA.sub_fold, idx_allAA_all, ...
			idx_matchedAA, covariates, all_cov, all_y, corr_mat, opt, metrics);
        yp_matchedAA(f) = yp;
        yt_matchedAA(f) = yt;
        for midx = 1:length(metrics)
			optimal_stats.(metrics{midx})(f,1) = pred_stats(midx);
        end
    end
    fprintf('\n')
	save(fullfile(model_dir, bhvr_nm{b}, ['final_result_matchedAA_' bhvr_nm{b} '.mat']), 'optimal_acc', 'optimal_stats', ...
		'yp_matchedAA', 'yt_matchedAA', 'y_matchedAA', 'y_matchedAA_resid', 'y_train_resid')
end
	
end