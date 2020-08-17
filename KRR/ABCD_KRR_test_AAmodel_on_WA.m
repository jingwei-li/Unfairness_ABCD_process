function ABCD_KRR_test_AAmodel_on_WA(csvname, model_dir, bhvr_ls, cfds_ls, full_subj_ls, ...
	full_FC, split_dir, split_fstem, AA_subdir, AAWA_stem)

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
	% load fold splits of AA and WA
	AA_fold = load(fullfile(split_dir, AA_subdir, [bhvr_nm{b} split_fstem '.mat']));
	AAWA_fold = load(fullfile(split_dir, [AAWA_stem split_fstem '_' bhvr_nm{b} '.mat']));
	allAA = CBIG_text2cell(fullfile(split_dir, AA_subdir, ['subj_' bhvr_nm{b} split_stem '.txt']));
	Nsplits = length(AA_fold.sub_fold);
	if(length(AAWA_fold.sub_fold) ~= Nsplits)
		error('Lengths of AA_fold and AAWA_fold are not equal.')
	end

	% load covariates of AA
	load(fullfile(model_dir, bhvr_nm{b}, ['confounds_' strjoin(cfds_nm, '_') '.mat']));
	if(size(covariates,2) ~= ncfds)
		error('Number of confounds in cfds_ls and model confounds .mat file not equal.')
	end

	% load optimal parameters
	opt = load(fullfile(model_dir, bhvr_nm{b}, ['final_result_' bhvr_nm{b} '.mat']))

	% collect raw behavioral scores of full set of subjects
	all_y = d.(bhvr_nm{b});
	all_y = all_y(idx_full_subj);

	[~, ~, idxAA_all] = intersect(allAA, all_subj, 'stable');
	yp_WA = cell(Nsplits, 1);  yt_WA = yp_WA;
	y_WA = yp_WA;  y_WA_resid = yp_WA;
	for f = 1:Nsplits
		fprintf('%d..', f)
		[~, ~, idx_WA] = intersect(AAWA_fold.sub_fold(f).selWA, all_subj, 'stable');

		[yp, yt, optimal_acc(f), pred_stats, y_WA{f}, y_WA_resid{f}] = ...
			ABCD_KRR_test_grp1Model_on_grp2_perfold( f, model_dir, bhvr_nm{b}, ...
			AA_fold.sub_fold, idxAA_all, idx_WA, covariates, all_cov, all_y, ...
			corr_mat, opt, metrics);
		yp_WA(f) = yp;
		yt_WA(f) = yt;
		for midx = 1:length(metrics)
			optimal_stats.(metrics{midx})(f,1) = pred_stats(midx);
		end
	end
	fprintf('\n')
	save(fullfile(model_dir, bhvr_nm{b}, ['final_result_randWA_' bhvr_nm{b} '.mat']), ...
		'optimal_acc', 'optimal_stats', 'yp_WA', 'yt_WA', 'y_WA', 'y_WA_resid')
end

end

