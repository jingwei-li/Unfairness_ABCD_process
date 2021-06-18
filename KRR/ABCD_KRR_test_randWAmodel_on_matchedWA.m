function ABCD_KRR_test_randWAmodel_on_matchedWA(csvname, model_dir, bhvr_ls, cfds_ls, cfds_X_ls, full_subj_ls, ...
	full_FC, split_dir, split_stem, WA_subdir)

% ABCD_KRR_test_randWAmodel_on_matchedWA(csvname, model_dir, bhvr_ls, cfds_ls, full_subj_ls, ...
%     full_FC, split_dir, split_stem, WA_subdir)
% 
% Check the performance of KRR models trained on randomly selected WA when they are tested on matched WA.
%
% Inputs:
% - csvname
%   Full path of the csv file containing all confounds and behaviors. It is generated by
%   `../preparation/ABCD_read_all_measures.m`. Default:
%   '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/phenotypes_pass_rs.txt'.
%
% - model_dir
%   Full path of the directory containing the kernel ridge regression models which were trained
%   on randomly selected WA. It is the output directory of `ABCD_KRR_reg_AgeSexMtIcvPEduc_y_randWA.sh`.
%
% - bhvr_ls
%   List of behavioral measures (full path). If you only want to test for one behavioral measure, 
%   this behavioral name can be directly passed in as a string. Default:
%   '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'.
%
% - cfds_ls
%   List of confounding variables which need to be regressed out from behavioral scores (full path).
%   Default: '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/confounds_list.txt'.
%
% - cfds_X_ls
%   List of confounding variables which need to be regressed out from RSFC (full path). Default: NONE.
%
% - full_subj_ls
%   List of subjects who have passed all quality controls and had all required phenotypes (full path).
%   It is generated by `../preparation/ABCD_read_all_measures.m`. Default:
%   '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt'.
%
% - full_FC
%   Full path of a .mat file containing a 3D matrix which is the resting-state functional 
%   connectivity of all the subjects in `full_subj_ls`. It is the output file of 
%   `../preparation/ABCD_check_RSFC_NaN.m`. Default:
%   '/home/jingweil/storage/MyProject/fairAI/ABCD_race/mat/RSFC/pass_rs_pass_pheno_5351.mat'.
%
% - split_dir
%
% - split_stem
%
% - WA_subdir
%   The split folds of the random WA should have been saved in a subfolder under `split_dir`. `WA_subdir` 
%   is the relative name of the subfolder.
%
% Author: Jingwei Li

%% default input arguments
[csvname, subj_hdr, d, bhvr_nm, nbhvr, cfds_nm, ncfds, cfds_X_nm, ncfds_X, full_subj_ls, all_subj, idx_full_subj, corr_mat] = ...
    ABCD_KRR_test_grp1Model_on_grp2_parse_args(csvname, bhvr_ls, cfds_ls, cfds_X_ls, full_subj_ls, full_FC);
    
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
	
	all_cov_y = CBIG_read_y_from_csv( {csvname}, subj_hdr, cfds_nm, cfds_types, full_subj_ls, 'NONE', ',' );
end

%% load confounds (to be regressed from RSFC)
fprintf('Loading confounds to be regressed from RSFC of all subjects ...\n')
if(strcmpi(cfds_X_ls, 'none'))
    all_cov_X = [];
else
    cfds_types = repmat({'continuous'}, 1, ncfds);
	sex_idx = strcmpi(cfds_nm, 'sex');
	if(any(sex_idx))
		cfds_types{sex_idx} = 'categorical';
	end
    all_cov_X = CBIG_read_y_from_csv( {csvname}, subj_hdr, cfds_X_nm, cfds_types, full_subj_ls, 'NONE', ',' );
end

metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
for b = 1:nbhvr
	fprintf('#%d behavior: %s\n', b, bhvr_nm{b});
	split_matched = load(fullfile(split_dir, ['sub_fold' split_stem '_' bhvr_nm{b} '.mat']));
	split_randWA = load(fullfile(split_dir, WA_subdir, [bhvr_nm{b} split_stem '.mat']));
	randWA = CBIG_text2cell(fullfile(split_dir, WA_subdir, ['subj_' bhvr_nm{b} split_stem '.txt']));
    Nsplits = length(split_randWA.sub_fold);
    if(length(split_matched.sub_fold) ~= Nsplits)
		error('Lengths of split_matched and split_allAA are not equal.')
    end
    
    % load covariates of the random WA
	load(fullfile(model_dir, bhvr_nm{b}, ['confounds_' strjoin(cfds_nm, '_') '.mat']));
	if(size(covariates, 2) ~= ncfds)
		error('Number of confounds in cfds_ls and model confounds .mat file not equal.')
    end
    
    % load optimal parameters
    opt = load(fullfile(model_dir, bhvr_nm{b}, ['final_result_' bhvr_nm{b} '.mat']));
    
    % collect raw behavioral scores of full set of subjects
	all_y = d.(bhvr_nm{b});
    all_y = all_y(idx_full_subj);
    
    [~, ~, idx_randWA_all] = intersect(randWA, all_subj, 'stable');
	yp_matchedWA = cell(Nsplits, 1);  yt_matchedWA = yp_matchedWA;
	y_matchedWA = yp_matchedWA;   y_matchedWA_resid = yp_matchedWA;   y_matchedAA_resid = yp_matchedWA;
	optimal_acc = nan(Nsplits, 1);
    fprintf('Fold = ')
    for f = 1:Nsplits
        fprintf('%d..', f)
		[~, ~, idx_matchedWA] = intersect(split_matched.sub_fold(f).selWA, all_subj, 'stable');
		[~, ~, idx_matchedAA] = intersect(split_matched.sub_fold(f).selAA, all_subj, 'stable');

		[yp, yt, optimal_acc(f), pred_stats, y_matchedWA{f}, y_matchedWA_resid{f}, y_train_resid{f}, y_matchedAA_resid{f}] = ...
			ABCD_KRR_test_grp1Model_on_grp2_perfold( f, model_dir, bhvr_nm{b}, split_randWA.sub_fold, idx_randWA_all, ...
			idx_matchedWA, covariates, all_cov_y, all_cov_X, all_y, corr_mat, opt, metrics, idx_matchedAA);

		ssr = sum((yt{1} - yp{1}).^2) ./ length(yt{1});
		sst = sum(([y_matchedAA_resid{f}; y_matchedWA_resid{f}] - mean(y_train_resid{f})).^2) ...
			./ length([y_matchedAA_resid{f}; y_matchedWA_resid{f}]);
		pred_stats(length(metrics) + 1) = bsxfun(@minus, 1, ssr./sst);
		fprintf('ssr = %f, sst = %f\n', ssr, sst)

        yp_matchedWA(f) = yp;
        yt_matchedWA(f) = yt;
        for midx = 1:length(metrics)
			optimal_stats.(metrics{midx})(f,1) = pred_stats(midx);
		end
		optimal_stats.pCOD(f,1) = pred_stats(length(metrics) + 1);
    end
    fprintf('\n')
	save(fullfile(model_dir, bhvr_nm{b}, ['final_result_matchedWA_' bhvr_nm{b} '.mat']), 'optimal_acc', 'optimal_stats', ...
		'yp_matchedWA', 'yt_matchedWA', 'y_matchedWA', 'y_matchedWA_resid', 'y_train_resid', 'y_matchedAA_resid')
end

end

