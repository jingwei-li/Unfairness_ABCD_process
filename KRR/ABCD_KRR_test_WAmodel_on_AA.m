function  ABCD_KRR_test_WAmodel_on_AA(csvname, model_dir, bhvr_ls, cfds_ls, full_subj_ls, ...
	full_FC, split_dir, split_fstem, WA_subdir, AAWA_subdir, outstem)

% ABCD_KRR_test_WAmodel_on_AA(csvname, model_dir, bhvr_ls, cfds_ls, full_subj_ls, ...
%	 full_FC, split_dir, split_fstem, WA_subdir, AAWA_subdir, outstem)
%
% Check the performance of KRR models trained on WA when they are tested on AA subjects.
%
% Inputs:
%   - csvname
%     Name of csv file containing all confounds and behaviors (full path).
%   - model_dir
%     Directory of models trained on WA (full path). It contains a subfolder for each behavior.
%   - bhvr_ls (optional)
%     List of behaviors (full path). Default:
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
%   - cfds_ls (optional)
%     List of confounding variable names (full path). Default:
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/confounds_list.txt'
%   - full_subj_ls (optional)
%     List of all subjects involved is project. Default:
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt'
%   - full_FC (optional)
%     Functional connectivity file corresponding to 'full_subj_ls'. Default:
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/mat/RSFC/pass_rs_pass_pheno_5351.mat'
%   - split_dir
%     Full path of the directory containing fold splits.
%   - split_fstem
%     Filename stem of fold splits.
%   - WA_subdir
%     Relative name of subdirectory containing fold splits of only WA subjects.
%   - AAWA_subdir
%     Relative name of subdirectory containing fold splits of the AA to be tested.
%   - outstem
%     Stem of output filename. The output filename will be 
%     fullfile(model_dir, <behavior>, ['final_result' outstem '_' <behavior> '.mat'])

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
	WA_fold = load(fullfile(split_dir, WA_subdir, [bhvr_nm{b} split_fstem '.mat']));
	AAWA_name = fullfile(split_dir, AAWA_subdir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']);
	if(~exist(AAWA_name, 'file'))
		AAWA_name = fullfile(split_dir, AAWA_subdir, [bhvr_nm{b} split_fstem '.mat']);
		if(~exist(AAWA_name, 'file'))
			error('Cannot find sub_fold file for WA based on split_dir = %s and AAWA_subdir = %s.', split_dir, AAWA_subdir)
		end
	end
	AAWA_fold = load(AAWA_name);
	WA = CBIG_text2cell(fullfile(split_dir, WA_subdir, ['subj_' bhvr_nm{b} split_fstem '.txt']));
	Nsplits = length(WA_fold.sub_fold);
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
    
    [~, ~, idxWA_all] = intersect(WA, all_subj, 'stable');
	yp_AA = cell(Nsplits, 1);  yt_AA = yp_AA;
	y_AA = yp_AA;  y_AA_resid = yp_AA;   y_matchedWA_resid = yp_AA;
	for f = 1:Nsplits
		fprintf('%d..', f)
		[~, ~, idx_AA] = intersect(AAWA_fold.sub_fold(f).selAA, all_subj, 'stable');
		[~, ~, idx_matchedWA] = intersect(AAWA_fold.sub_fold(f).selWA, all_subj, 'stable');

		[yp, yt, optimal_acc(f), pred_stats, y_AA{f}, y_AA_resid{f}, y_train_resid{f}, y_matchedWA_resid{f}] = ...
			ABCD_KRR_test_grp1Model_on_grp2_perfold( f, model_dir, bhvr_nm{b}, ...
			WA_fold.sub_fold, idxWA_all, idx_AA, covariates, all_cov, all_y, ...
			corr_mat, opt, metrics, idx_matchedWA);

		ssr = sum((yt{1} - yp{1}).^2) ./ length(yt{1});
		sst = sum(([y_AA_resid{f}; y_matchedWA_resid{f}] - mean(y_train_resid{f})).^2) ...
			./ length([y_AA_resid{f}; y_matchedWA_resid{f}]);
		pred_stats(length(metrics) + 1) = bsxfun(@minus, 1, ssr./sst);
		fprintf('ssr = %f, sst = %f\n', ssr, sst)

		yp_AA(f) = yp;
		yt_AA(f) = yt;
		for midx = 1:length(metrics)
			optimal_stats.(metrics{midx})(f,1) = pred_stats(midx);
		end
		optimal_stats.pCOD(f,1) = pred_stats(length(metrics) + 1);
	end
	fprintf('\n')
	save(fullfile(model_dir, bhvr_nm{b}, ['final_result' outstem '_' bhvr_nm{b} '.mat']), ...
		'optimal_acc', 'optimal_stats', 'yp_AA', 'yt_AA', 'y_AA', 'y_AA_resid', 'y_train_resid', 'y_matchedWA_resid')
end
    
end