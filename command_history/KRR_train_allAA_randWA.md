# Train KRR on all AA and random WA (see `KRR_allAA_vs_randWA.md`)

## Create `sub_fold` structure for selected subjects

```matlab
proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
ABCD_create_subfold_allAA_randWA(...
    fullfile(proj_dir, 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt'), ...
    fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
    fullfile(proj_dir, 'mat', 'matchANDsplit', '20200719'), ...
    '_pass_rs_pass_pheno')
```

## ------ KRR: regress age, sex, FD, DVARS, ICV, parental education from behaviors ------

### Train KRR

```bash
../KRR/ABCD_KRR_reg_AgeSexMtIcvPEduc_y_allAA_randWA.sh
```

### Test the all-AA + random-WA trained model with matched WA subjects

```bash
../KRR/ABCD_KRR_test_allAArandWAmodel_on_matchedWA.sh
```

### Compare accuracy of matched WA between models trained on whole population and trained on allAA+random WA

1. metric: predictive COD

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('All-AA & random-WA', 'WA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA_randWA'), 'predictive_COD', ...
        'Compare matched WA accuracy: whole-population trained VS all-AA + random-WA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_randWA_trained', 'KRR', '20200721'), 'pCOD_matchedWA_WholePopModel_vs_allAArandWAmodel')
    ```

2. metric: Pearson's correlation

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('All-AA & random-WA', 'WA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA_randWA'), 'corr', ...
        'Compare matched WA accuracy: whole-population trained VS all-AA + random-WA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_randWA_trained', 'KRR', '20200721'), 'corr_matchedWA_WholePopModel_vs_allAArandWAmodel')
    ```

### Test the all-AA + random-WA trained model with matched AA subjects

```bash
../KRR/ABCD_KRR_test_allAArandWAmodel_on_matchedAA.sh
```

### Compare accuracy of matched AA between models trained on whole population and trained on allAA+random WA

1. metric: predictive COD

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('All-AA & random-WA', 'AA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA_randWA'), 'predictive_COD', ...
        'Compare matched AA accuracy: whole-population trained VS all-AA + random-WA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_randWA_trained', 'KRR', '20200721'), 'pCOD_matchedAA_WholePopModel_vs_allAArandWAmodel')
    ```

2. metric: Pearson's correlation

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('All-AA & random-WA', 'AA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA_randWA'), 'corr', ...
        'Compare matched AA accuracy: whole-population trained VS all-AA + random-WA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_randWA_trained', 'KRR', '20200721'), 'corr_matchedAA_WholePopModel_vs_allAArandWAmodel')
    ```

### Using model trained on all AA and random WA, compare accuracy between matched AA and matched WA

1. metric: predictive COD

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA_randWA'), ...
        'matchedAA', 'matchedWA', [], [], 'predictive_COD', ...
        'Matched AA vs WA (model trained on all AA + random WA)', ...
        fullfile(proj_dir, 'figures', 'allAA_randWA_trained', 'KRR', '20200721'), ...
        'pCOD_allAArandWAmodel_test_matchedAAvsWA')
    ```

2. metric: Pearson's correlation

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA_randWA'), ...
        'matchedAA', 'matchedWA', [], [], 'corr', ...
        'Matched AA vs WA (model trained on all AA + random WA)', ...
        fullfile(proj_dir, 'figures', 'allAA_randWA_trained', 'KRR', '20200721'), ...
        'corr_allAArandWAmodel_test_matchedAAvsWA')
    ```

3. metric: COD

    ```matlab
    proj_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA_randWA'), ...
        'matchedAA', 'matchedWA', [], [], 'COD', ...
        'Matched AA vs WA (model trained on all AA + random WA)', ...
        fullfile(proj_dir, 'figures', 'allAA_randWA_trained', 'KRR', '20200721'), ...
        'COD_allAArandWAmodel_test_matchedAAvsWA')
    ```
