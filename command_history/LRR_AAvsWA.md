## ------ LRR: regress age, sex, FD, DVARS, ICV, parental education from both behaviors and RSFC ------

### Run LRR (bash)

```bash
../LRR/ABCD_LRR_reg_AgeSexMtIcvPEduc_fr_y_fc.sh
```

### Plot correlation and predictive COD accuracy of all behaviors on all test subjects

```matlab
```

### Permutation test of predictability

#### Step 1. generate permutaions

```matlab
ABCD_LRR_predictable_behavior_step1('/home/jingweil/storage/MyProject/fairAI/ABCD_race/models/LRR/20210629/reg_AgeSexMtIcvPEduc_fr_y_FC', 1000)
```

#### Step2. repeat LRR with permutated y for each behavior

```bash
```