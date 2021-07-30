## Create RSFC without censoring

1. Generate lists used for `CBIG_preproc_FCmetrics.m`

```matlab
ABCD_FC_no_censor_genlists('/home/jingweil/storage/MyProject/fairAI/ABCD_race/mat/RSFC/no_censor')
```

2. Compute RSFC for each subject

```bash
../FC_no_censor/ABCD_FC_no_censor.sh
```

3. Concatenate across subjects, Fisher's z-transformation

```matlab
FC_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/mat/RSFC/no_censor';
ABCD_FC_no_censor_concat(FC_dir, fullfile(FC_dir, 'individuals'))
```

## Run kernel ridge regression

```bash
../KRR/ABCD_KRR_reg_AgeSexMtIcvPEduc_from_y_FC.sh -FC_file /home/jingweil/storage/MyProject/fairAI/ABCD_race/mat/RSFC/no_censor/RSFC_5351_no_censor_z.mat -outdir /home/jingweil/storage/MyProject/fairAI/ABCD_race/models/KRR/20200721/FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC
```

### Plot correlation of all behaviors on all test subjects

```matlab
proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
ABCD_KRR_violin_acc_allsub(fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
   'corr', fullfile(proj_dir, 'figures', 'AAvsWA', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
   'corr_allsubj_violin')
```

### Permutation test of predictability

1. metric: predictive COD

   ```bash
   ssh headnode
   proj_dir=/home/jingweil/storage/MyProject/fairAI/ABCD_race
   cmd="matlab -nodesktop -nodisplay -nojvm -r \" cd('$proj_dir/scripts/Unfairness_ABCD_process');\
      ABCD_addpath; ABCD_KRR_predictable_behavior(\
      fullfile('$proj_dir', 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), \
      fullfile('$proj_dir', 'mat', 'matchANDsplit', '20200719'), \
      '_pass_rs_pass_pheno', 1000, 'predictive_COD', \
      fullfile('$proj_dir', 'mat', 'predictability', 'KRR', 'pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'));  exit; \" "
   echo $cmd
   work_dir=$proj_dir/models/KRR/20200721/FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC/logs
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
      fullfile('$proj_dir', 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), \
      fullfile('$proj_dir', 'mat', 'matchANDsplit', '20200719'), \
      '_pass_rs_pass_pheno', 1000, 'corr', \
      fullfile('$proj_dir', 'mat', 'predictability', 'KRR', 'corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'));  exit; \" "
   echo $cmd
   work_dir=$proj_dir/models/KRR/20200721/FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC/logs
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
model_dir = fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC');
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
      fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt'), ...
      fullfile(proj_dir, 'mat', 'matchANDsplit', '20200719'), ...
      '_pass_rs_pass_pheno', 120, ...
      fullfile(proj_dir, 'mat', 'predictability', 'KRR', 'pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'R_thres0.15_12behaviors.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'))
   ```

2. metric: Pearson's correlation

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_KRR_corr_AAvsWA(...
      fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt'), ...
      fullfile(proj_dir, 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt'), ...
      fullfile(proj_dir, 'mat', 'matchANDsplit', '20200719'), ...
      '_pass_rs_pass_pheno', 120, ...
      fullfile(proj_dir, 'mat', 'predictability', 'KRR', 'corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'R_thres0.15_12behaviors.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'))
   ```

### Permutation test of accuracy difference

1. metric: predictive COD

   Only predictable behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable.txt'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'sig_pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat') )
   ```

   All behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      [], ...
      'predictive_COD', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'sig_pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat') )
   ```

2. metric: Pearson's correlation

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'corr_predictable.txt'), ...
      'corr', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'sig_corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat') )
   ```

   All behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_PermTest_AAvsWA( ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      [], ...
      'corr', ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'sig_corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat') )
   ```

### Violin plot

1. metric: predictive COD

   Only plot predictable behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_AAvsWA(...
      fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', ...
      'pCOD_predictable.txt'),  fullfile(proj_dir, 'models', 'KRR', '20200721', ...
      'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'pCOD_predictable_colloquial.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'sig_pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), 'pCOD')
   ```

   Plot all behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_AAvsWA([], [], ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'sig_pCOD_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      'predictive_COD', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), 'pCOD_allbehaviors')
   ```

2. metric: Pearson's correlation
   
   Only plot predictable behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_AAvsWA(...
      fullfile(proj_dir, 'models', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', ...
      'corr_predictable.txt'),  fullfile(proj_dir, 'models', 'KRR', '20200721', ...
      'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC', 'lists', 'corr_predictable_colloquial.txt'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'sig_corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC_predictable.mat'), ...
      'corr', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), 'corr')
   ```

   Plot all behaviors:

   ```matlab
   proj_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race';
   ABCD_violin_AAvsWA([], [], ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      fullfile(proj_dir, 'mat', 'AAvsWA', 'KRR', 'sig_corr_pass_rs_pass_pheno_FCnocensor_reg_AgeSexMtIcvPeduc_fr_y_FC.mat'), ...
      'corr', ...
      fullfile(proj_dir, 'figures', 'AAvsWA', 'KRR', '20200721', 'FCnocensor_reg_AgeSexMtIcvPEduc_fr_y_FC'), 'corr_allbehaviors')
   ```
