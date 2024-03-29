# Train KRR on all AA

## Create `sub_fold` structure for all AA subjects

```matlab
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
ABCD_create_subfold_allAA(...
    fullfile(proj_dir, 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt'), ...
    fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
    fullfile(proj_dir, 'mat', 'matchANDsplit', '20200719'), ...
    '_pass_rs_pass_pheno')
```

## ------ KRR: regress age, sex, FD, DVARS, ICV, parental education from behaviors and RSFC ------

### Run KRR only on all AA

```bash
cd ../KRR
./ABCD_KRR_reg_AgeSexMtIcvPEduc_y_FC_allAA.sh
```

### Test the all-AA model with matched AA subjects

```bash
proj_dir=/home/jingweil/storage/MyProject/fairAI/ABCD_race
../KRR/ABCD_KRR_test_AllAAmodel_on_matchedAA.sh -model_dir \
  $proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_FC_allAA -cfds_X_ls \
  /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/confounds_list.txt
```

### Test the all-AA model with matched WA

```bash
proj_dir=/home/jingweil/storage/MyProject/fairAI/ABCD_race
../KRR/ABCD_KRR_test_AllAAmodel_on_matchedWA.sh -model_dir \
  $proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_FC_allAA -cfds_X_ls \
  /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/confounds_list.txt
```

### Using model trained on AA, compare accuracy between matched AA and matched WA

1. metric: predictive COD

    Permutation test

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_FC_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', [], 'predictive_COD', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'pCOD_allAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA.mat'))

    % behaviors that are predictable in the whole-population model
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_FC_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable.txt'), ...
        'predictive_COD', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'pCOD_allAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA_predictable_WholepopModel.mat'))
    ```

    Plotting

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_FC_allAA'), ...
        'matchedAA', 'matchedWA', [], [], 'predictive_COD', ...
        'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'pCOD_AllAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'pCOD_allAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA.mat'))

    % plot only the behaviors that are predictable in the whole-population model
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_FC_allAA'), ...
        'matchedAA', 'matchedWA', ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable.txt'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable_colloquial.txt'), ...
        'predictive_COD', ...
        'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'pCOD_AllAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA_predictable_WholepopModel', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'pCOD_allAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA_predictable_WholepopModel.mat'))
    ```

2. metric: Pearson's correlation

    Permutation test

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_FC_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', [], 'corr', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'corr_allAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA.mat'))

    % behaviors that were predictable in the whole-population model
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_FC_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'corr_predictable.txt'), ...
        'corr', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'corr_allAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA_predictable_WholepopModel.mat'))
    ```

    Plotting

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_FC_allAA'), ...
        'matchedAA', 'matchedWA', [], [], 'corr', ...
        'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'corr_AllAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'corr_allAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA.mat'))

    % Plot only the behaviors that are predictable in the whole-population model
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_FC_allAA'), ...
        'matchedAA', 'matchedWA', ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'corr_predictable.txt'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'corr_predictable_colloquial.txt'), ...
        'corr', 'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'corr_AllAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA_predictable_WholepopModel', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'corr_allAAmodel_reg_AgeSexMtIcvPEduc_y_FC_test_matchedAAvsWA_predictable_WholepopModel.mat'))
    ```

### Plot chord diagrams of model-learned brain-behavioral association

```matlab
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
bhvr_ls = fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt');
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
for b = 1:nbhvr
    ABCD_chord_learnedBBA(fullfile(proj_dir, 'mat', 'interpretation', 'KRR', '20200719', ...
        'reg_AgeSexMtIcvPEduc_fr_y_FC_allAA', ['learned_BWAS_' bhvr_nm{b} '.mat']), ...
        fullfile(proj_dir, 'figures', 'interpretation', 'KRR', '20200721', ...
        'reg_AgeSexMtIcvPEduc_fr_y_FC_allAA', ['chord_' bhvr_nm{b}]));
end
```

## ------ KRR: regress age, sex, FD, DVARS, ICV, parental education from behaviors ------

### Test the all-AA model with matched AA subjects

```bash
../KRR/ABCD_KRR_test_AllAAmodel_on_matchedAA.sh
```

### Compare accuracy of matched AA between models trained on whole population and trained on all AA

1. metric: predictive COD

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('AA', 'AA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'predictive_COD', ...
        'Compare matched AA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'pCOD_matchedAA_WholePopModel_vs_AllAAmodel')
    ```

2. metric: Pearson's correlation

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('AA', 'AA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'corr', ...
        'Compare matched AA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'corr_matchedAA_WholePopModel_vs_AllAAmodel')
    ```

3. metric: MSE

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('AA', 'AA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'MSE_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'MSE', ...
        'Compare matched AA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'MSE_matchedAA_WholePopModel_vs_AllAAmodel')
    ```

### Test the all-AA model with matched WA

```bash
../KRR/ABCD_KRR_test_AllAAmodel_on_matchedWA.sh
```

### Compare accuracy of matched WA between models trained on whole population and trained on all AA

1. metric: predictive COD

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('AA', 'WA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'predictive_COD', ...
        'Compare matched WA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'pCOD_matchedWA_WholePopModel_vs_AllAAmodel')
    ```

2. metric: Pearson's correlation

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('AA', 'WA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'corr', ...
        'Compare matched WA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'corr_matchedWA_WholePopModel_vs_AllAAmodel')
    ```

3. metric: MSE
   
    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('AA', 'WA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'MSE_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'MSE', ...
        'Compare matched WA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'MSE_matchedWA_WholePopModel_vs_AllAAmodel')
    ```

4. metric: normalized MSE

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whiskerYA_trainAll_vs_trainXA('AA', 'WA', [], [], ...
        fullfile(proj_dir, 'mat', 'AAvsWA', 'normMSE_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y.mat'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), 'MSE_norm', ...
        'Compare matched WA accuracy: whole-population trained VS all-AA trained', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), 'normMSE_matchedWA_WholePopModel_vs_AllAAmodel')
    ```

### Using model trained on AA, compare accuracy between matched AA and matched WA

1. metric: predictive COD

    Permutation test

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', [], 'predictive_COD', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'pCOD_allAAmodel_test_matchedAAvsWA.mat'))

    % behaviors that are predictable in the whole-population model
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y', 'lists', 'pCOD_predictable.txt'), ...
        'predictive_COD', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'pCOD_allAAmodel_test_matchedAAvsWA_predictable_WholepopModel.mat'))
    ```

    Plotting

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedWA', [], [], 'predictive_COD', ...
        'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'pCOD_AllAAmodel_test_matchedAAvsWA', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'pCOD_allAAmodel_test_matchedAAvsWA.mat'))

    % plot only the behaviors that are predictable in the whole-population model
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedWA', ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y', 'lists', 'pCOD_predictable.txt'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y', 'lists', 'pCOD_predictable_colloquial.txt'), ...
        'predictive_COD', ...
        'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'pCOD_AllAAmodel_test_matchedAAvsWA_predictable_WholepopModel', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'pCOD_allAAmodel_test_matchedAAvsWA_predictable_WholepopModel.mat'))
    ```

2. metric: Pearson's correlation

    Permutation test

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', [], 'corr', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'corr_allAAmodel_test_matchedAAvsWA.mat'))

    % behaviors that were predictable in the whole-population model
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y', 'lists', 'corr_predictable.txt'), ...
        'corr', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'corr_allAAmodel_test_matchedAAvsWA_predictable_WholepopModel.mat'))
    ```

    Plotting

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedWA', [], [], 'corr', ...
        'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'corr_AllAAmodel_test_matchedAAvsWA', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'corr_allAAmodel_test_matchedAAvsWA.mat'))

    % Plot only the behaviors that are predictable in the whole-population model
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedWA', ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y', 'lists', 'corr_predictable.txt'), ...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_fr_y', 'lists', 'corr_predictable_colloquial.txt'), ...
        'corr', 'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'corr_AllAAmodel_test_matchedAAvsWA_predictable_WholepopModel', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'corr_allAAmodel_test_matchedAAvsWA_predictable_WholepopModel.mat'))
    ```

3. metric: MSE

    Permutation test

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_PermTest_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedAA', 'matchedWA', 'WA', [], 'MSE', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'MSE_allAAmodel_test_matchedAAvsWA.mat'))
    ```

    Plotting

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedWA', [], [], 'MSE', ...
        'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'MSE_AllAAmodel_test_matchedAAvsWA', ...
        fullfile(proj_dir, 'mat', 'AAvsWA_train_subpop', 'MSE_allAAmodel_test_matchedAAvsWA.mat'))
    ```

4. metric: COD

    ```matlab
    proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
    ABCD_KRR_whisker_trainXA_testYAvsZA(...
        fullfile(proj_dir, 'models', 'KRR', '20200721', 'reg_AgeSexMtIcvPEduc_y_allAA'), ...
        'matchedAA', 'matchedWA', [], [], 'COD', ...
        'Matched AA vs WA (model trained on all AA)', ...
        fullfile(proj_dir, 'figures', 'allAA_trained', 'KRR', '20200721'), ...
        'COD_AllAAmodel_test_matchedAAvsWA')
    ```
