function ABCD_FC_subgroup(full_subj_ls, sgrp_subj_ls, full_FC, out_FC)

%

full_subj = CBIG_text2cell(full_subj_ls);
sgrp_subj = CBIG_text2cell(sgrp_subj_ls);

[~, ~, idx] = intersect(sgrp_subj, full_subj, 'stable');
load(full_FC)
corr_mat = corr_mat(:,:,idx);
fprintf('FC of subgroup extracted.\n')

outdir = fileparts(out_FC);
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(out_FC, 'corr_mat', '-v7.3')
fprintf('Saved.\n')

end