
#

```matlab
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
ABCD_cmp_betas_dist(fullfile(proj_dir, 'figures', 'cmp_beta_orig_vs_resample'), ...
    fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
    fullfile(proj_dir, 'same_size_AAWA', 'betas'))
```