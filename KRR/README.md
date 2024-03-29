## Case 1. Training on whole population

### Run kernel ridge regression

The top-level scripts are `ABCD_KRR_reg_AgeSexMtIcvPEduc_from_y.sh` and `ABCD_KRR_reg_AgeSexMtIcvPEduc_from_FC.sh`. These two scripts handle the cases where confounding variables (age, gender, FD, DVARS, intracranial volume, parental education) are regressed from behavioral scores and from resting-state functional connectivity respectively. Exemplar usage of `ABCD_KRR_reg_AgeSexMtIcvPEduc_from_y.sh` (bash):

```bash
./ABCD_KRR_reg_AgeSexMtIcvPEduc_from_y.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -subfold_dir /.../matchANDsplit/ \
    -subfold_stem _pass_rs_pass_pheno -outdir /.../reg_AgeSexMtIcvPEduc_fr_y/
```

`<path_to_lists>/phenotypes_pass_rs.txt`, `<path_to_lists>/behavior_list.txt`, `<path_to_lists>/subjects_pass_rs_pass_pheno.txt` are output files of `../preparation/ABCD_read_all_measures.m`. `<path_to_lists>/confounds_list.txt` is a list specified by the user containing the confounding variables to be regressed out. `<path_to_RSFC>/pass_rs_pass_pheno_5351.mat` is the output of `../preparation/ABCD_check_RSFC_NaN.m`. `/.../matchANDsplit/` is the output directory of `../match_split/ABCD_match_and_split.m`. `_pass_rs_pass_pheno` should be consistent with the `outstem` argument passed into `../match_split/ABCD_match_and_split.m`.

### Check which behavioral measures showed prediction accuracy significantly above chance

Apply multi-level block permutation test to check whether the out-of-sample prediction accuracy across all subjects (all ethnic/racial groups, not only AA and WA) was significantly above chance. Exemplar usage (matlab):

```matlab
ABCD_KRR_predictable_behavior('/path/to/KRR/output/', '/path/to/split/folds/', ...
    '_pass_rs_pass_pheno', 1000, 'predictive_COD', '/output/predictability.mat', ...
    '/path/to/lists/behavior_list.txt', '/path/to/lists/colloquial_names_list.txt', ...
    '/path/to/lists/subjects_pass_rs_pass_pheno.txt', )
```

`'/path/to/KRR/output/'` is the output directory of `ABCD_KRR_reg_AgeSexMtIcvPEduc_from_y.sh` or `ABCD_KRR_reg_AgeSexMtIcvPEduc_from_FC.sh`. `'/path/to/split/folds/'` is the output directory of `../match_split/ABCD_match_and_split.m`. `'_pass_rs_pass_pheno'` is a string attached to the output files, specified by the user. `1000` is the number of multi-level block permutations. `'predictive_COD'` is the accuracy metric. `'/path/to/lists/behavior_list.txt'` and `'/path/to/lists/colloquial_names_list.txt'` are the list of behavioral names and their colloquial names generated by `../preparation/ABCD_read_all_measures.m`. `'/path/to/lists/subjects_pass_rs_pass_pheno.txt'` is the list of subjects who passed all quality controls and had all required phenotypes, generated by `../preparation/ABCD_read_all_measures.m`. `'/path/to/lists/phenotypes_pass_rs.txt'` is the csv file generated by `../preparation/ABCD_read_all_measures.m`.

## Case 2. Training on all AA and random WA with the same sample size (i.e. 50% AA, 50% WA)

### Run KRR on the subpopulation

```bash
./ABCD_KRR_reg_AgeSexMtIcvPEduc_y_allAA_randWA.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -full_subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -full_FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -split_dir /.../matchANDsplit/train_allAA_randWA \
    -split_fstem _pass_rs_pass_pheno -outdir /.../reg_AgeSexMtIcvPEduc_y_allAA_randWA
```

`/.../matchANDsplit/train_allAA_randWA` is the directory containing split folds of all AA + random WA. It is the output directory of `../match_split/ABCD_create_subfold_allAA_randWA.m`. `/.../reg_AgeSexMtIcvPEduc_y_allAA_randWA` is the output directory to store the KRR models trained on all AA + random WA.

### Test the subpopulation-trained model on matched AA

Calculate predicted behavioral scores and prediction accuracy of matched AA for the KRR models trained only on all AA + random WA.

```bash
./ABCD_KRR_test_allAArandWAmodel_on_matchedAA.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -full_subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -full_FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -split_dir /.../matchANDsplit/ \
    -split_fstem _pass_rs_pass_pheno -train_folds_subdir train_allAA_randWA
```

`/.../matchANDsplit/` is the directory containing split folds of the whole population. It is the output directory of `../match_split/ABCD_match_and_split.m`. `train_allAA_randWA` is the name of a subdirectory under `/.../matchANDsplit/` which contains split folds of all AA + random WA. This subdirectory is the output directory of `../match_split/ABCD_create_subfold_allAA_randWA.m`.

### Test the subopopulation-trained model on matched WA

Calculate predicted behavioral scores and prediction accuracy of matched WA for the KRR models trained only on all AA + random WA.

Use `ABCD_KRR_test_allAArandWAmodel_on_matchedWA.sh` to calculate out-of-sample prediction accuracy of matched WA when KRR models were trained on all AA + random WA. Input arguments of `ABCD_KRR_test_allAArandWAmodel_on_matchedWA.sh` are the same with that of `ABCD_KRR_test_allAArandWAmodel_on_matchedAA.sh`.

## Case 3: training solely on AA

### Run KRR on the subpopulation

```bash
./ABCD_KRR_reg_AgeSexMtIcvPEduc_y_allAA.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -full_subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -full_FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -split_dir /.../matchANDsplit/train_allAA \
    -split_fstem _pass_rs_pass_pheno -outdir /.../reg_AgeSexMtIcvPEduc_y_allAA
```

`/.../matchANDsplit/train_allAA` is the directory containing split folds of all AA. It is the output directory of `../match_split/ABCD_create_subfold_allAA.m`. `/.../reg_AgeSexMtIcvPEduc_y_allAA` is the output directory to store the KRR models trained on all AA.

### Test the subpopulation-trained model on matched AA

Calculate predicted behavioral scores and prediction accuracy of matched AA for the KRR models trained only on all AA.

```bash
./ABCD_KRR_test_AllAAmodel_on_matchedAA.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -full_subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -full_FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -split_dir /.../matchANDsplit/ \
    -split_fstem _pass_rs_pass_pheno -AA_subfolder train_allAA
```

`train_allAA` is the name of a subdirectory under `/.../matchANDsplit/` which contains split folds of all AA. This subdirectory is the output directory of `../match_split/ABCD_create_subfold_allAA.m`.

### Test the subopopulation-trained model on matched WA

Calculate predicted behavioral scores and prediction accuracy of matched WA for the KRR models trained only on all AA.

```bash
./ABCD_KRR_test_AllAAmodel_on_matchedWA.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -full_subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -full_FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -split_dir /.../matchANDsplit/ \
    -split_fstem _pass_rs_pass_pheno -AA_subdir train_allAA
```

## Case 4: training solely on WA

### Run KRR on the subpopulation

```bash
./ABCD_KRR_reg_AgeSexMtIcvPEduc_y_randWA.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -full_subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -full_FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -split_dir /.../matchANDsplit/train_randWA \
    -split_fstem _pass_rs_pass_pheno -outdir /.../reg_AgeSexMtIcvPEduc_y_randWA
```

`/.../matchANDsplit/train_randWA` is the directory containing split folds of all AA + random WA. It is the output directory of `../match_split/ABCD_create_subfold_randWA.m`. `/.../reg_AgeSexMtIcvPEduc_y_randWA` is the output directory to store the KRR models trained on random WA.

### Test the subpopulation-trained model on matched AA

Calculate predicted behavioral scores and prediction accuracy of matched AA for the KRR models trained only on random WA.

```bash
./ABCD_KRR_test_randWAmodel_on_matchedAA.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -full_subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -full_FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -split_dir /.../matchANDsplit/ \
    -split_fstem _pass_rs_pass_pheno -WA_subfolder train_randWA
```

`train_randWA` is the name of a subdirectory under `/.../matchANDsplit/` which contains split folds of random WA. This subdirectory is the output directory of `../match_split/ABCD_create_subfold_randWA.m`.

### Test the subopopulation-trained model on matched WA

Calculate predicted behavioral scores and prediction accuracy of matched WA for the KRR models trained only on random WA.

```bash
./ABCD_KRR_test_randWAmodel_on_matchedWA.sh -csvname <path_to_lists>/phenotypes_pass_rs.txt \
    -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls <path_to_lists>/confounds_list.txt \
    -full_subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt -full_FC_file \
    <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -split_dir /.../matchANDsplit/ \
    -split_fstem _pass_rs_pass_pheno -WA_subdir train_randWA
```
