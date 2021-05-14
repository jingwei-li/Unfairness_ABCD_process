#!/bin/sh
#
# Save out predicted behavioral scores of training subjects.
# Author: Jingwei Li

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir="/home/jingweil/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y"
has_subdir=0
split_dir="$proj_dir/mat/matchANDsplit/20200719"
split_fstem="_pass_rs_pass_pheno"

main() {
    log_dir="$model_dir/logs"
    jname="training_prediction"

    cmd="matlab -nodesktop -nojvm -nodisplay -r \" addpath('$DIR'); \
        ABCD_KRR_training_prediction('$model_dir', $has_subdir, '$split_dir', \
        '$split_fstem'); exit;  \" "

    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 02:00:00 -mem 12G \
        -name $jname -joberr $log_dir/$jname.err -jobout $log_dir/$jname.out
}

#############################
# Function usage
#############################
usage() { GREEN='\033[0;32m'; NOCOLOR='\033[0m'; echo -e "
NAME:
    ABCD_KRR_training_prediction_reg_AgeSexMtIcvPEduc_fr_y.sh

DESCRIPTION:
    Save out predicted behavioral scores of training subjects.

ARGUMENTS:
    -h                         : Display help information.
    -model_dir   <model_dir>   : Full path of the directory containing kernel ridge regression
                                 models trained on whole population. It is the output directory
                                 of ../KRR/ABCD_KRR_reg_AgeSexMtIcvPEduc_from_*.sh. Default:
                                 ${GREEN}$model_dir${NOCOLOR}
    -split_dir   <split_dir>   : Full path of the directory containing the split folds of whole
                                 population. It is the output directory of 
                                 ../match_split/ABCD_match_and_split.m. Default:
                                 ${GREEN}$split_dir${NOCOLOR}
    -split_fstem <split_fstem> : A string shared across all filenames of split folds. It should 
                                 be the same string as passed into 
                                 ../match_split/ABCD_match_and_split.m. Default: 
                                 ${GREEN}$split_fstem${NOCOLOR}
" 1>&2; exit 1; }

##########################################
# Parse Arguments 
##########################################

while [[ $# -gt 0 ]]; do
    flag=$1; shift;

    case $flag in
		-h) usage; exit ;;
        -model_dir) model_dir=$1; shift ;;
        -split_dir) split_dir=$1; shift ;;
        -split_fstem) split_fstem=$1; shift ;;
        *) echo "Unknown flag: $flag"; usage; 1>&2; exit 1 
    esac
done

main
