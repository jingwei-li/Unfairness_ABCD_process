function [csvname, subj_hdr, d, bhvr_nm, nbhvr, cfds_nm, ncfds, full_subj_ls, all_subj, idx_full_subj, corr_mat] = ...
	ABCD_KRR_test_grp1Model_on_grp2_parse_args(csvname, bhvr_ls, cfds_ls, full_subj_ls, full_FC)

% [csvname, d, bhvr_nm, nbhvr, cfds_nm, ncfds, all_subj, full_subj_ls, corr_mat] = ...
%	 ABCD_KRR_test_grp1Model_on_grp2_parse_args(csvname, bhvr_ls, cfds_ls, full_subj_ls, full_FC)
%
% Parse input arguments for the functions which test the models trained on grp1 using the subjects
% in grp2.
%
% Inputs (all paths should be full paths):
%   - csvname: Name of csv file containing all behavioral and confonding variables.
%   - bhvr_ls: Behaviors list. Default:
%         '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
%   - cfds_ls: Confounds list. Default:
%         '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/confounds_list.txt'
%   - full_subj_ls: List of all subjects involved is project. Default:
%         '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt'
%   - full_FC: Functional connectivity file corresponding to 'full_subj_ls'. Default:
%         '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/mat/RSFC/pass_rs_pass_pheno_5351.mat'
%
% Outputs:
%   - csvname: return default value if input 'csvname' is empty.
%   - d: table read from 'csvname'.
%   - bhvr_nm: behavioral names.
%   - nbhvr: number of behaviors.
%   - cfds_nm: confounding variable names.
%   - ncfds: number of confounds.
%   - full_subj_ls: return default value if input 'full_subj_ls' is empty.
%   - all_subj: cell of subject IDs read from 'full_subj_ls'.
%   - corr_mat: #features x #subjects functional connectivity, read from 'full_FC'.
%               #features is (#ROI - 1) x #ROI / 2; #subjects should be the length
%               of 'all_subj'.
%

proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
ls_dir = fullfile(proj_dir, 'scripts', 'lists');
if(~exist('csvname', 'var') || isempty(csvname))
	csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
	bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
if(ischar(bhvr_ls))
	if(exist(bhvr_ls, 'file'))
		[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
	else
		bhvr_nm = {bhvr_ls}; nbhvr = 1;
	end
elseif(iscell(bhvr_ls))
	bhvr_nm = bhvr_ls;  nbhvr = length(bhvr_nm);
end

if(~exist('cfds_ls', 'var') || isempty(cfds_ls))
	cfds_ls = fullfile(ls_dir, 'confounds_list.txt');
end
[cfds_nm, ncfds] = CBIG_text2cell(cfds_ls);
if(~exist('full_subj_ls', 'var') || isempty(full_subj_ls))
	full_subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
all_subj = CBIG_text2cell(full_subj_ls);

if(~exist('full_FC', 'var') || isempty(full_FC))
	full_FC = fullfile(proj_dir, 'mat', 'RSFC', 'pass_rs_pass_pheno_5351.mat');
end
load(full_FC);
tril_ind = find(tril(ones(size(corr_mat,1), size(corr_mat,2)), -1) == 1);
corr_mat = reshape(corr_mat, size(corr_mat,1)*size(corr_mat,2), size(corr_mat,3));
corr_mat = corr_mat(tril_ind,:);
if(size(corr_mat) ~= length(all_subj))
	error('full_subj_ls and full_FC do not correspond to each other.')
end

subj_hdr = 'subjectkey';
d = readtable(csvname);
[~, ~, idx_full_subj] = intersect(all_subj, d.(subj_hdr), 'stable');

end