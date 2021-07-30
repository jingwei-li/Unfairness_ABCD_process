function ABCD_resample_betas(outdir, bhvr_nm, subj_ls, sub_fold_mat, csvname, cfds_ls, cfds_X_ls, orig_subj_ls, RSFC_file)

% ABCD_resample_betas()
%
% 

proj_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts');

if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(proj_dir, 'lists', 'phenotypes_pass_rs.txt');
end

if(~exist('cfds_ls', 'var') || isempty(cfds_ls))
    cfds_ls = fullfile(proj_dir, 'lists', 'confounds_list.txt');
end
if(strcmpi(cfds_ls, 'none'))
    cfds_nm = {'none'};
else
    [cfds_nm, Ncfds] = CBIG_text2cell(cfds_ls);
end

if(~exist('cfds_X_ls', 'var') || isempty(cfds_X_ls))
    cfds_X_ls = 'NONE';
end
if(strcmpi(cfds_X_ls, 'none'))
    cfds_X_nm = {'none'};
else
    [cfds_X_nm, Ncfds_X] = CBIG_text2cell(cfds_X_ls);
end

if(~exist('orig_subj_ls', 'var') || isempty(orig_subj_ls))
    orig_subj_ls = fullfile(proj_dir, 'lists', 'subjects_pass_rs_pass_pheno.txt');
end

if(~exist('RSFC_file', 'var') || isempty(RSFC_file))
    RSFC_file = fullfile(proj_dir, 'mat', 'RSFC', 'pass_rs_pass_pheno_5351.mat');
end

[subjects, nsub] = CBIG_text2cell(subj_ls);
orig_subj = CBIG_text2cell(orig_subj_ls);
[~, idx] = ismember(subjects, orig_subj);
subj_hdr = 'subjectkey';

%% grab behavior, write to a .mat file
y_file = fullfile(outdir, ['y_' bhvr_nm '.mat']);
if(~exist(y_file, 'file'))
    y = CBIG_read_y_from_csv( {csvname}, subj_hdr, {bhvr_nm}, {'continuous'}, subj_ls, y_file, ',' );
else
    load(y_file)
end

%% grab covariates which need to be regressed from behaviors, save to a .mat file
cfds_file = fullfile(outdir, ['confounds_' strjoin(cfds_nm, '_'), '.mat']);
if(~exist(cfds_file, 'file'))
    if(strcmpi(cfds_ls, 'none'))
        fprintf('No regressor to be regressed from behaviors.\n')
        covariates = 'NONE';
    else
        cfds_types = repmat({'continuous'}, 1, Ncfds);
        sex_idx = strcmpi(cfds_nm, 'sex');
        if(any(sex_idx))
            cfds_types{sex_idx} = 'categorical';
        end
        
        covariates = CBIG_read_y_from_csv( {csvname}, subj_hdr, cfds_nm, cfds_types, subj_ls, 'NONE', ',' );
    end
    save(cfds_file, 'covariates')
else
    load(cfds_file)
end

%% grab covariates which need to be regressed from RSFC, save to a .mat file
cfds_X_file = fullfile(outdir, ['confounds_X_' strjoin(cfds_X_nm, '_'), '.mat']);
if(~exist(cfds_X_file, 'file'))
    if(strcmpi(cfds_X_ls, 'none'))
        fprintf('No regressor to be regressed from RSFC.\n')
        cov_X = [];
    else
        cfds_types = repmat({'continuous'}, 1, Ncfds_X);
        sex_idx = strcmpi(cfds_nm, 'sex');
        if(any(sex_idx))
            cfds_types{sex_idx} = 'categorical';
        end

        cov_X = CBIG_read_y_from_csv( {csvname}, subj_hdr, cfds_X_nm, cfds_types, subj_ls, 'NONE', ',' );
    end
    save(cfds_X_file, 'cov_X')
end

%% Calculate betas for regression from behavioral scores
load(sub_fold_mat)
CBIG_crossvalid_regress_covariates_from_y( y, covariates, sub_fold, outdir, bhvr_nm);

%% Calculate betas for regression from RSFC
Nfolds = length(sub_fold);
load(RSFC_file)
corr_mat = corr_mat(:,:,idx);
tril_ind = find(tril(ones(size(corr_mat,1)), -1) == 1);
corr_mat = reshape(corr_mat, size(corr_mat,1)*size(corr_mat,2), size(corr_mat,3));
corr_mat = corr_mat(tril_ind, :);
for f = 1:Nfolds
    train_ind = sub_fold(f).fold_index == 0;
    test_ind = sub_fold(f).fold_index == 1;
    
    % Regress covariates from features
    curr_outdir = fullfile(outdir, 'FSM_innerloop', ['fold_' num2str(f)]);
    if(~exist(curr_outdir, 'dir')); mkdir(curr_outdir); end
    if(exist('cov_X', 'var') && ~isempty(cov_X))
        fprintf('Fold %d, regress covariates from features\n', f)
        [feature_resid(train_ind, :), beta] = CBIG_regress_X_from_y_train(...
            corr_mat(:, train_ind)', cov_X(train_ind, :));

        save(fullfile(curr_outdir, 'beta.mat'), 'beta')
    end
end

end