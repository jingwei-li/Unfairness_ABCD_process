function ABCD_regress_cfds_from_FC(FC_in, csvname, cfds_ls, subj_ls, FC_out)

% ABCD_regress_cfds_from_FC(FC_in, csvname, cfds_ls, subj_ls, FC_out)
%
% Input:
% - FC_in
%   Input RSFC filename (full path).
%
% - csvname (optional)
%   Filename of CSV which contains the confound values (full path).
%   Default:
%   '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/phenotypes_pass_rs.txt'
% 
% - cfds_ls (optional)
%   List of confounds (full path). Default:
%   '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/confounds_list.txt'
% 
% - subj_ls (optional)
%   List of subjects (full path). Default:
%   '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt'
%   
% - FC_out
%   Output RSFC filename (full path).
%


%% default arguments
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end

if(~exist('cfds_ls', 'var') || isempty(cfds_ls))
    cfds_ls = fullfile(ls_dir, 'confounds_list.txt');
end
[cfds_nm, Ncfds] = CBIG_text2cell(cfds_ls);

if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
subj_hdr = 'subjectkey';


%% load
cfds_types = repmat({'continuous'}, 1, Ncfds);
sex_idx = strcmpi(cfds_nm, 'sex');
if(any(sex_idx))
    cfds_types{sex_idx} = 'categorical';
end
covariates = CBIG_read_y_from_csv( {csvname}, subj_hdr, cfds_nm, cfds_types, subj_ls, 'NONE', ',' );

load(FC_in)
if(size(corr_mat,3) ~= size(covariates,1))
    error('Number of subjects in FC and subj_ls are different!')
end

%% regress
% reshape FC from 3D to 2D
s = size(corr_mat);
corr_mat = reshape(corr_mat, s(1)*s(2), s(3));

% regression: regressors are demeaned
corr_mat = CBIG_glm_regress_matrix(corr_mat', ...
    bsxfun(@minus, covariates, mean(covariates,1)), -1, []);

% reshape back
corr_mat = reshape(corr_mat', s);

% save
outdir = fileparts(FC_out);
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(FC_out, 'corr_mat', '-v7.3')

end

