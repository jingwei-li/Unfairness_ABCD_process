## Save out predicted behavioral scores of training subjects

This step facilitates the next step, since model-learned brain-behavioral associations are computated as the covariance between RSFC and predicted behavioral scores in the **training** set. Use:

```bash
./ABCD_KRR_training_prediction_reg_AgeSexMtIcvPEduc_fr_y.sh -model_dir \
    /path/to/KRR/models/  -split_dir /.../matchANDsplit/ -split_fstem \
    _pass_rs_pass_pheno
```

## Calculate model-learned brain-behavioral association patterns

```bash
./ABCD_KRR_learned_BWAS_reg_AgeSexMtIcvPEduc_fr_y.sh -model_dir \
    /path/to/KRR/models/  -split_dir /.../matchANDsplit/ -split_fstem \
    _pass_rs_pass_pheno -outdir /output/mat/dir/ -figdir /output/figure/dir/
```

## Calculate true brain-behavioral association patterns of matched AA and WA separately

True brain-behavioral associations are computed as the covariance between RSFC and true behavioral scores in the test set within each ethnic/racial group. Use:

```bash
./ABCD_KRR_true_BWAS_AAvsWA_reg_AgeSexMtIcvPEduc_fr_y.sh -model_dir \
    /path/to/KRR/models/  -split_dir /.../matchANDsplit/ -split_fstem \
    _pass_rs_pass_pheno -outdir /output/mat/dir/ -figdir /output/figure/dir/
```
