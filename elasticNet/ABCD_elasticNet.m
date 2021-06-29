function ABCD_elasticNet(csvname, bhvr_nm, cfds_ls, cfds_X_ls, subj_ls, subfold_f, FC_file, N_inner_folds, outdir, outstem)

% ABCD_elasticNet(csvname, bhvr_nm, cfds_ls, cfds_X_ls, subj_ls, subfold_f, FC_file, N_inner_folds, outdir, outstem)
%
% Wrapper function to apply CBIG elastic net package in the ABCD dataset.
%
% Inputs:
% - csvname
%   Full path of the csv file containing all phenotypical values which are necessary for this study.
%   It is generated by `../preparation/ABCD_read_all_measures.m`.
%   Default: '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/phenotypes_pass_rs.txt'
%
% - bhvr_nm
%   Behavioral measure to be predicted.
%
% - cfds_ls
%   List of confounding variables which need to be regressed from behavioral measures (full path).
%   Default: '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/confounds_list.txt'
%
% - cfds_X_ls 
%   List of confounding variables which need to be regressed from RSFC (full path). Pass in 'NONE' if
%   nothing needs to be regressed from RSFC.
%
% - subj_ls
%   List of subjects who have passed all quality controls and had all required phenotypes (full path).
%   It is generated by `../preparation/ABCD_read_all_measures.m`.
%   Default: '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt'
%
% - subfold_f
%   Filename of split folds, full path. It is the output of `../match_split/ABCD_match_and_split.m`.
%
% - FC_file
%   Functional connectivity filename, full path.
%
% - N_inner_folds
%   The number of hyperparameter-selection (inner-loop) cross-validation folds. 
%
% - outdir
%   Output directory, full path.
% 
% - outstem
%   A discrimitive string to be attached to output filenames.
%
% Author: Jingwei Li

ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');

if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end

if(~exist('cfds_ls', 'var') || isempty(cfds_ls))
    cfds_ls = fullfile(ls_dir, 'confounds_list.txt');
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

if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_ls);
subj_hdr = 'subjectkey';

if(ischar(N_inner_folds))
    N_inner_folds = str2double(N_inner_folds);
end

%% grab behavior, write to a .mat file
y_file = fullfile(outdir, ['y_' bhvr_nm '.mat']);
if(~exist(y_file, 'file'))
    CBIG_read_y_from_csv( {csvname}, subj_hdr, {bhvr_nm}, {'continuous'}, subj_ls, y_file, ',' );
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

%% Call elastic net workflow
params = load(subfold_f, 'sub_fold');
load(FC_file); params.feature_mat = corr_mat; clear corr_mat
load(cfds_file); params.covariates = covariates; clear covariates
load(cfds_X_file); params.cov_X = cov_X; clear cov_X
load(y_file); params.y = y; clear y
params.outdir = outdir;
params.outstem = bhvr_nm; params.num_inner_folds = N_inner_folds; params.split_name = [];

CBIG_run_Elasticnet_workflow(params)
    
end