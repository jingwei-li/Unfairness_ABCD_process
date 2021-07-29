function ABCD_FC_truncate_concat( outdir, indiv_dir, subj_ls )

% 

if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt';
end

[subjects, nsub] = CBIG_text2cell(subj_ls);
prefix = 'RSFC_5351_truncate';

corr_mat = [];
for s = 1:nsub
    fname = fullfile(indiv_dir, subjects{s}, [prefix '_all2all.mat']);
    FC = load(fname);
    corr_mat = cat(3, corr_mat, FC.corr_mat);
end
save(fullfile(outdir, [prefix, '.mat']), 'corr_mat', '-v7.3')
corr_mat = CBIG_StableAtanh(corr_mat);
save(fullfile(outdir, [prefix, '_z.mat']), 'corr_mat', '-v7.3')

end
