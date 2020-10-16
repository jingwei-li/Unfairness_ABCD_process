#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo $DIR

proj_dir="$HOME/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y"
has_subdir=0
split_dir="$proj_dir/mat/matchANDsplit/20200719"
split_fstem="_pass_rs_pass_pheno"
outmat="$proj_dir/mat/interpretation/KRR/20200719/reg_AgeSexMtIcvPEduc_fr_y/learned_BWAS.mat"
figdir="$proj_dir/figures/AAvsWA/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y/interpretations/learned_BWAS"

log_dir="$proj_dir/mat/interpretation/KRR/20200719/reg_AgeSexMtIcvPEduc_fr_y"
mkdir -p $log_dir
jname="ABCD_learned_BWAS"
cmd="matlab -nodesktop -nodisplay -r \" addpath('$DIR'); \
    ABCD_KRR_learned_BWAS('$model_dir', $has_subdir, '$split_dir', '$split_fstem', \
    '$outmat', '$figdir'); exit; \" "

echo $cmd
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 10:00:00 -mem 30G \
    -name $jname -joberr $log_dir/$jname.err -jobout $log_dir/$jname.out