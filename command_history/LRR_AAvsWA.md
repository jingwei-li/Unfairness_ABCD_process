## ------ LRR: regress age, sex, FD, DVARS, ICV, parental education from both behaviors and RSFC ------

### Run LRR (bash)

```bash
../LRR/ABCD_LRR_reg_AgeSexMtIcvPEduc_fr_y_fc.sh
```

### Plot correlation and predictive COD accuracy of all behaviors on all test subjects

```matlab
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
ABCD_LRR_violin_acc_allsub(fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
   'corr', fullfile(proj_dir, 'figures', 'AAvsWA', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
   'corr_allsubj_violin')
```

### Permutation test of predictability

#### Step 1. generate permutaions

```matlab
ABCD_LRR_predictable_behavior_step1('/home/jingweil/storage/MyProject/fairAI/ABCD_race/models/LRR/20210629/reg_AgeSexMtIcvPEduc_fr_y_FC', 1000)
```

#### Step2. repeat LRR with permuted y for each behavior

```bash
../LRR/ABCD_LRR_predictable_behavior_step2.sh -Nperm 1000
```

#### Step 3: compute p values

1. Accuracy metric: predictive COD

```matlab
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
ABCD_LRR_predictable_behavior_step3(...
    fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), 'predictive_COD', ...
    fullfile(proj_dir, 'mat', 'predictability', 'LRR', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'), ...
    fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
    fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt'), 1000)
```

2. Accuracy metric: Pearson's correlation

```matlab
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
ABCD_LRR_predictable_behavior_step3(...
    fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), 'corr', ...
    fullfile(proj_dir, 'mat', 'predictability', 'LRR', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'), ...
    fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
    fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt'), 1000)
```

### Compute accuracy metric per race group; get predictable behaviors

Get behavioral measures with >0.15 correlation accuracy across all test subjects:

```matlab
clear
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
bhvr_ls = fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt');
colloq_ls = fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt');
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
colloq_nm = CBIG_text2cell(colloq_ls);
model_dir = fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC');
opt = zeros(nbhvr, 1);
for b = 1:nbhvr
    load(fullfile(model_dir, bhvr_nm{b}, 'results', 'optimal_acc', [bhvr_nm{b} '.mat']))
    Nfolds = length(optimal_statistics);
    curr_corr = zeros(Nfolds,1);
    for f = 1:Nfolds
        curr_corr(f) = optimal_statistics{f}.corr;
    end
    opt(b) = mean(curr_corr);
end
idx = find(opt > 0.1);
bhvr_nm = bhvr_nm(idx);
colloq_nm = colloq_nm(idx);
mkdir(fullfile(model_dir, 'lists'))
CBIG_cell2text(bhvr_nm, fullfile(model_dir, 'lists', ['R_thres0.1_' num2str(length(idx)) 'behaviors.txt']))
CBIG_cell2text(colloq_nm, fullfile(model_dir, 'lists', ['R_thres0.1_' num2str(length(idx)) 'colloquial.txt']))
```

1. metric: predictive COD

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_LRR_pCOD_AAvsWA(...
      fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt'), ...
      fullfile(proj_dir, 'mat', 'matchANDsplit', '20200719'), ...
      '_pass_rs_pass_pheno', 120, ...
      fullfile(proj_dir, 'mat', 'predictability', 'LRR', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'R_thres0.1_12behaviors.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'))
   ```

2. metric: Pearson's correlation

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_LRR_corr_AAvsWA(...
      fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt'), ...
      fullfile(proj_dir, 'mat', 'matchANDsplit', '20200719'), ...
      '_pass_rs_pass_pheno', 120, ...
      fullfile(proj_dir, 'mat', 'predictability', 'LRR', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'R_thres0.1_12behaviors.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'))
   ```

### Permutation test of accuracy difference

1. metric: predictive COD

   Only predictable behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable.txt'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'sig_pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat') )
   ```

   All behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      [], ...
      'predictive_COD', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'sig_pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat') )
   ```

2. metric: Pearson's correlation

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'corr_predictable.txt'), ...
      'corr', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'sig_corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat') )
   ```

   All behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      [], ...
      'corr', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'sig_corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat') )
   ```

### Violin plot

1. metric: predictive COD

   Only plot predictable behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_withnull_AAvsWA(...
      fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', ...
      'pCOD_predictable.txt'),  fullfile(proj_dir, 'models', 'LRR', '20210629', ...
      'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable_colloquial.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'sig_pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), 'pCOD')
   ```

   Plot all behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_withnull_AAvsWA([], [], ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'sig_pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), 'pCOD_allbehaviors')
   ```

2. metric: Pearson's correlation
   
   Only plot predictable behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_withnull_AAvsWA(...
      fullfile(proj_dir, 'models', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', ...
      'corr_predictable.txt'),  fullfile(proj_dir, 'models', 'LRR', '20210629', ...
      'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'corr_predictable_colloquial.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'sig_corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      'corr', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), 'corr')
   ```

   Plot all behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_withnull_AAvsWA([], [], ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', 'sig_corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      'corr', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), 'corr_allbehaviors')
   ```

#### Plot "predicted score  - original score"

All 36 behaviors

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';

   % Are AA-WA differences in terms of "predicted - orginal score" significant?
   ABCD_PermTest_predVStrue_AAvsWA(fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', ...
       'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), [], ...
       fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', ...
       'sig_predVStrue_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'))

   % plot
   ABCD_KRR_violin_predVStrue(...
       fullfile(proj_dir, 'figures', 'AAvsWA', 'LRR', '20210629', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
       'predVStrue_allbehaviors', fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', ...
       'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
       fullfile(proj_dir, 'mat', 'AAvsWA', 'LRR', ...
       'sig_predVStrue_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), [])
   ```

### Collect optimal parameters

```matlab
ABCD_LRR_collect_opt_param('/home/jingweil/storage/MyProject/fairAI/ABCD_race/models/LRR/20210629/reg_AgeSexMtIcvPEduc_fr_y_FC', [])
```