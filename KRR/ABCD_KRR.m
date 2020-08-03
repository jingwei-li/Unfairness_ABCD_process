function ABCD_KRR(csvname, bhvr_nm, cfds_ls, subj_ls, subfold_f, FC_file, N_inner_folds, outdir, outstem)

ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';

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

%% grab covariates, save to a .mat file
cfds_file = fullfile(outdir, ['confounds_' strjoin(cfds_nm, '_'), '.mat']);
if(~exist(cfds_file, 'file'))
    if(strcmpi(cfds_ls, 'none'))
        fprintf('No regressor to be regressed from behavios.\n')
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

%%


CBIG_KRR_workflow_LITE( [], 0, subfold_f, y_file, ...
    cfds_file, FC_file, N_inner_folds, outdir, outstem )

end

