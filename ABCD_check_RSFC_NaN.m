function pass_subj = ABCD_check_RSFC_NaN(fmri_dir, subj_list)

if(~exist('fmri_dir', 'var') || isempty(fmri_dir))
    fmri_dir = '/mnt/eql/yeo13/data/ABCD/rs_GSR';
end

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt';
end
[subjects, nsub] = CBIG_text2cell(subj_list);

hasnan = logical(zeros(nsub,1));
for s = 1:nsub
    RSFC_file = fullfile(fmri_dir, subjects{s}, 'FC_metrics', 'Pearson_r', [subjects{s} ...
        '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08_fs6_sm6_all2all.mat']);
    load(RSFC_file)
    hasnan(s) = any(isnan(corr_mat(:)));
    clear corr_mat
end

pass_subj = subjects(~hasnan);

end

