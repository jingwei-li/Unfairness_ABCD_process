# Train KRR on all AA

## ------ KRR: regress age, sex, FD, DVARS, ICV, parental education from behaviors ------

### Compare accuracy of matched AA between models trained on whole population and trained on all AA

1. metric: predictive COD

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_whiskerYA_trainAll_vs_trainXA('AA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'predictive_COD', ...
        'Compare matched AA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'pCOD_matchedAA_WholePopModel_vs_AllAAmodel')
    ```

2. metric: Pearson's correlation

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_whiskerYA_trainAll_vs_trainXA('AA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'corr', ...
        'Compare matched AA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'corr_matchedAA_WholePopModel_vs_AllAAmodel')
    ```

### Compare accuracy of matched WA between models trained on whole population and trained on all AA

1. metric: predictive COD

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_whiskerYA_trainAll_vs_trainXA('WA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'predictive_COD', ...
        'Compare matched WA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'pCOD_matchedWA_WholePopModel_vs_AllAAmodel')
    ```

2. metric: Pearson's correlation

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_whiskerYA_trainAll_vs_trainXA('WA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'corr', ...
        'Compare matched WA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'corr_matchedWA_WholePopModel_vs_AllAAmodel')
    ```

### Using model trained on AA, compare accuracy between matched AA and matched WA

