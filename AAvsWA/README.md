## Compare out-of-sample prediction accuracy between matched AA and WA when knernel ridge regression models were trained on whole population

### Calculate accuracy of matched AA and WA

Predictive COD, Pearson's correlation, and prediction error metrics can be calculated using `ABCD_KRR_pCOD_AAvsWA.m`, `ABCD_KRR_corr_AAvsWA.m` and `ABCD_KRR_error_AAvsWA.m` respectively. Taking `ABCD_KRR_pCOD_AAvsWA.m` as an example:

```matlab
ABCD_KRR_pCOD_AAvsWA('/whole/population/KRR/dir/', '/path/to/lists/behavior_list.txt', ...
    '/path/to/lists/colloquial_list.txt', '/path/to/lists/subjects_pass_rs_pass_pheno.txt', ...
    '/.../matchANDsplit/', '_pass_rs_pass_pheno', 120, '/path/to/predictable/behaviors.mat', ...
    '/behavioral/list/with/high-enough/overall/acc.txt', '/path/to/output.mat', 0)
```

`'/whole/population/KRR/dir/'` is the output directory of `../KRR/ABCD_KRR_reg_AgeSexMtIcvPEduc_from_*.sh`. `'/path/to/lists/behavior_list.txt'`, `'/path/to/lists/colloquial_list.txt'`, `'/path/to/lists/subjects_pass_rs_pass_pheno.txt'` are the output files of `../preparation/ABCD_read_all_measures.m`. `'/.../matchANDsplit/'` is the output directory of `../match_split/ABCD_match_and_split.m`. `'/path/to/predictable/behaviors.mat'` is the output file of `ABCD_KRR_predictable_behavior.m`. `'/behavioral/list/with/high-enough/overall/acc.txt'` is a list of behavioral measures which achieved a non-trivial prediction accuracy across all test subjects (including all ethnic/racial groups).

### Significant difference in prediction accuracy

```matlab
ABCD_PermTest_AAvsWA( '/accuracy/of/matchedAAWA.mat', '/path/to/lists/behavior_list.txt', ...
    'predictive_COD', '/output/significancy.mat' )
```

`'/accuracy/of/matchedAAWA.mat'` is the output file of `ABCD_KRR_pCOD_AAvsWA.m`.

## Compare out-of-sample prediction accuracy between all AA and random WA when knernel ridge regression models were trained on whole population

### Calculate accuracy of all AA and random WA

Predictive CPD and Pearson's correlation can be calculated using `ABCD_KRR_pCOD_allAA_vs_randWA.m` and `ABCD_KRR_corr_allAA_vs_randWA.m` respectively. Taking `ABCD_KRR_pCOD_allAA_vs_randWA.m` as an example:

```matlab
ABCD_KRR_pCOD_allAA_vs_randWA('/whole/population/KRR/dir/', '/path/to/lists/behavior_list.txt', ...
    '/path/to/lists/subjects_pass_rs_pass_pheno.txt', '/.../matchANDsplit/', ..
    '_pass_rs_pass_pheno', 120, '/path/to/output.mat')
```

### Significant difference in prediction accuracy

```matlab
ABCD_PermTest_AAvsWA( '/accuracy/of/allAA_vs_randomWA.mat', '/path/to/lists/behavior_list.txt', ...
    'predictive_COD', '/output/significancy.mat' )
```

`'/accuracy/of/allAA_vs_randomWA.mat'` is the output file of `ABCD_KRR_pCOD_AAvsWA.m`.
