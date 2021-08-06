# KRR: Only behavioral scores are matched between AA and WA

## Preparation

### Split & match

```matlab
ABCD_match_and_split(...
   '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/phenotypes_pass_rs.txt', ...
   '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt', ...
   'race', 'site', 'family_id', ...
   'none', ...
   '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt', 100, 2.45, ...
   '/home/jingweil/storage/MyProject/fairAI/ABCD_race/mat/matchANDsplit/20210802', ...
   '_pass_rs_pass_pheno')
```

## ------ KRR: regress age, sex, FD, DVARS, ICV, parental education from both behaviors and RSFC ------

### Run KRR (bash)

```bash
../KRR/ABCD_KRR_reg_AgeSexMtIcvPEduc_from_y_FC.sh -subfold_dir \
/home/jingweil/storage/MyProject/fairAI/ABCD_race/mat/matchANDsplit/20210802 \
-outdir /home/jingweil/storage/MyProject/fairAI/ABCD_race/models/KRR/20210802/reg_AgeSexMtIcvPEduc_fr_y_FC
```

### Plot correlation and predictive COD accuracy of all behaviors on all test subjects

```matlab
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
ABCD_KRR_violin_acc_allsub(fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
   'corr', fullfile(proj_dir, 'figures', 'AAvsWA', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
   'corr_allsubj_violin')
```

### Permutation test of predictability

1. metric: predictive COD

   ```bash
   ssh headnode
   proj_dir=/home/jingweil/storage/MyProject/fairAI/ABCD_race
   cmd="matlab -nodesktop -nodisplay -nojvm -r \" cd('$proj_dir/scripts/Unfairness_ABCD_process');\
      ABCD_addpath; ABCD_KRR_predictable_behavior(\
      fullfile('$proj_dir', 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), \
      fullfile('$proj_dir', 'mat', 'matchANDsplit', '20210802'), \
      '_pass_rs_pass_pheno', 1000, 'predictive_COD', \
      fullfile('$proj_dir', 'mat', 'predictability', 'KRR', '20210802', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'));  exit; \" "
   work_dir=$proj_dir/models/KRR/20210802/reg_AgeSexMtIcvPEduc_fr_y_FC/logs
   jname=perm_pCOD_predictability
   $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 30:00:00 -mem 12G \
      -name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
   ```

2. metric: Pearson's correlation
   
   ```bash
   ssh headnode
   proj_dir=/home/jingweil/storage/MyProject/fairAI/ABCD_race
   cmd="matlab -nodesktop -nodisplay -nojvm -r \" cd('$proj_dir/scripts/Unfairness_ABCD_process'); \
      ABCD_addpath; ABCD_KRR_predictable_behavior(\
      fullfile('$proj_dir', 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), \
      fullfile('$proj_dir', 'mat', 'matchANDsplit', '20210802'), \
      '_pass_rs_pass_pheno', 1000, 'corr', \
      fullfile('$proj_dir', 'mat', 'predictability', 'KRR', '20210802', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'));  exit; \" "
   work_dir=$proj_dir/models/KRR/20210802/reg_AgeSexMtIcvPEduc_fr_y_FC/logs
   jname=perm_corr_predictability
   $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 30:00:00 -mem 12G \
      -name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
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
model_dir = fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC');
opt = zeros(nbhvr, 1);
for b = 1:nbhvr
   load(fullfile(model_dir, bhvr_nm{b}, ['final_result_' bhvr_nm{b} '.mat']))
   opt(b) = mean(optimal_stats.corr);
end
idx = find(opt > 0.15);
bhvr_nm = bhvr_nm(idx);
colloq_nm = colloq_nm(idx);
mkdir(fullfile(model_dir, 'lists'))
CBIG_cell2text(bhvr_nm, fullfile(model_dir, 'lists', ['R_thres0.15_' num2str(length(idx)) 'behaviors.txt']))
CBIG_cell2text(colloq_nm, fullfile(model_dir, 'lists', ['R_thres0.15_' num2str(length(idx)) 'colloquial.txt']))
```

1. metric: predictive COD

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_KRR_pCOD_AAvsWA(...
      fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt'), ...
      fullfile(proj_dir, 'mat', 'matchANDsplit', '20210802'), ...
      '_pass_rs_pass_pheno', 120, ...
      fullfile(proj_dir, 'mat', 'predictability', 'KRR', '20210802', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'R_thres0.15_12behaviors.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'))
   ```

2. metric: Pearson's correlation

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_KRR_corr_AAvsWA(...
      fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt'), ...
      fullfile(proj_dir, 'mat', 'matchANDsplit', '20210802'), ...
      '_pass_rs_pass_pheno', 120, ...
      fullfile(proj_dir, 'mat', 'predictability', 'KRR', '20210802', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'R_thres0.15_12behaviors.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'))
   ```

### Permutation test of accuracy difference

1. metric: predictive COD

   Only predictable behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable.txt'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'sig_pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat') )
   ```

   All behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      [], ...
      'predictive_COD', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'sig_pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat') )
   ```

2. metric: Pearson's correlation

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'corr_predictable.txt'), ...
      'corr', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'sig_corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat') )
   ```

   All behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      [], ...
      'corr', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'sig_corr_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat') )
   ```

### Violin plot

#### Plot accuracy difference

1. metric: predictive COD

   Only plot predictable behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_AAvsWA(...
      fullfile(proj_dir, 'models', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', ...
      'pCOD_predictable.txt'),  fullfile(proj_dir, 'models', 'KRR', '20210802', ...
      'reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable_colloquial.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'sig_pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), 'pCOD')
   ```

   Plot all behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_AAvsWA([], [], ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', '20210802', 'sig_pCOD_pass_rs_pass_pheno_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'KRR', '20210802', 'reg_AgeSexMtIcvPEduc_fr_y_FC'), 'pCOD_allbehaviors')
   ```
