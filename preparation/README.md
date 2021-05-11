## Read (behavioral, demographic, morphological) phenotypes from ABCD csv files. Censor out the subjects with none of the required measures. Combine these measures into one csv file.

The user only needs to run the top-level script `ABCD_read_all_measures.m`. This script call other scriot such as `ABCD_read_gender.m`, `ABCD_read_LittleMan.m`, etc. Usage example (matlab):

```matlab
ABCD_read_all_measures('/preprocessed/rs-fMRI/dir/', '/sMRI/recon_all/dir/', ...
    '/list/of/subjects/passed/rs-fMRI/QC.txt', ...
    '/output/dir/', '<a shared string attached on all output files>')
```

## Sanity check

### Collected on Philips scanners?

Data of subjects which were collected on Philips scanners should be excluded, given the known issue: https://github.com/ABCD-STUDY/fMRI-cleanup. 

Use

```matlab
[subj_philips, isPhilips] = ABCD_get_subj_Philips(subj_list);
```

to check if any subject in `subj_list` was collected on Philips scanners.

### NaN in RSFC?

Use 

```matlab
pass_subj = ABCD_check_RSFC_NaN('/preprocessed/rs-fMRI/dir/', subj_list, '/concatenated/RSFC.mat');
```

to check if all subjects in `subj_list` have NaN values in the resting-state functional connectivity. `subj_list` is one of the output files of `ABCD_read_all_measures.m`. `'/preprocessed/rs-fMRI/dir/'` and `'/concatenated/RSFC.mat'` should be adjusted based on your data paths.

### Dependency of missing data on ethnicities/races?

Use 

```matlab
[nan_per_eth, nonnan_per_eth] = ABCD_missing_pheno_per_race(pheno_csv, bhvr_ls)
```

to count the number of empty and non-empty entries of each behavioral measure in `bhvr_ls`. Both `pheno_csv` and `bhvr_ls` are the output files of `ABCD_read_all_measures.m`.
