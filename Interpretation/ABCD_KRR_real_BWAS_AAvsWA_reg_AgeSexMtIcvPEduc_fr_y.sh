#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir="$HOME/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y"
has_subdir=0
split_dir="$proj_dir/mat/matchANDsplit/20200719"
split_fstem="_pass_rs_pass_pheno"
outmat="$proj_dir/mat/interpretation/KRR/20200719/reg_AgeSexMtIcvPEduc_fr_y/real_BWAS_testAAvsWA.mat"
figdir="$proj_dir/figures/AAvsWA/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y/interpretations/real_BWAS_testAAvsWA"

log_dir="$proj_dir/mat/interpretation/KRR/20200719/reg_AgeSexMtIcvPEduc_fr_y"
mkdir -p $log_dir
jname="ABCD_real_BWAS_testAAvsWA"

cmd="matlab -nodekstop -nodisplay -r \" addpath $DIR; \
    ABCD_real_BWAS_AAvsWA('$model_dir', $has_subdir, '$split_dir', '$split_fstem', '$outmat', '$fig_dir'); exit; \" "

echo $cmd
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 10:00:00 -mem 30G \
    -name $jname -joberr $log_dir/$jname.err -jobout $log_dir/$jname.out