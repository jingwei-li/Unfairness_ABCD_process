#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir="/home/jingweil/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y"
has_subdir=0
split_dir="$proj_dir/mat/matchANDsplit/20200719"
split_fstem="_pass_rs_pass_pheno"

log_dir="$model_dir/logs"
jname="training_prediction"

cmd="matlab -nodesktop -nojvm -nodisplay -r \" addpath('$DIR'); \
    ABCD_KRR_training_prediction('$model_dir', $has_subdir, '$split_dir', \
    '$split_fstem'); exit;  \" "

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 02:00:00 -mem 12G \
    -name $jname -joberr $log_dir/$jname.err -jobout $log_dir/$jname.out

