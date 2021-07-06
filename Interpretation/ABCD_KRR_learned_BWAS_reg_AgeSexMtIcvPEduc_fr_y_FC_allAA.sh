#!/bin/sh
#
# Calculate model-learned brain-behavioral association patterns
# Author: Jingwei Li

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo $DIR

proj_dir="$HOME/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_FC_allAA"
has_subdir=1
split_dir="$proj_dir/mat/matchANDsplit/20200719/train_allAA"
split_fstem="_pass_rs_pass_pheno"
outdir="$proj_dir/mat/interpretation/KRR/20200719/reg_AgeSexMtIcvPEduc_fr_y_FC_allAA"
figdir="$proj_dir/figures/AAvsWA/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y_FC_allAA/interpretations/learned_BWAS"

bhvr_ls="$proj_dir/scripts/lists/behavior_list.txt"
colloq_ls="$proj_dir/scripts/lists/colloquial_list.txt"
bhvr_nm=($(cat $bhvr_ls))

main() {
    log_dir=$outdir/logs
    mkdir -p $log_dir

    i=0
    while IFS= read -r colloq_nm; do
        echo $colloq_nm
        curr_bhvr=${bhvr_nm[$i]}
        if [ "$curr_bhvr" == "" ]; then continue; fi
        curr_FC="$model_dir/FC/$curr_bhvr.mat"
        outmat="$outdir/learned_BWAS_${curr_bhvr}.mat"

        jname="ABCD_learned_BWAS_$curr_bhvr"
        cmd="matlab -nodesktop -nodisplay -r \" addpath('$DIR'); \
            ABCD_KRR_learned_BWAS('$model_dir', $has_subdir, '$split_dir', '$split_fstem', \
            '$outmat', '$figdir', '$curr_FC', '$curr_bhvr', '$colloq_nm'); exit; \" "

        echo $cmd
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 01:00:00 -mem 20G \
            -name $jname -joberr $log_dir/$jname.err -jobout $log_dir/$jname.out
        i=$((i + 1))
    done < $colloq_ls
}

#############################
# Function usage
#############################
usage() { GREEN='\033[0;32m'; NOCOLOR='\033[0m'; echo -e "
NAME:
    ABCD_KRR_learned_BWAS_reg_AgeSexMtIcvPEduc_fr_y_FC.sh

DESCRIPTION:
    Calculate model-learned brain-behavioral association patterns

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
    -outdir      <outdir>      : Output directory of the main .mat file and log files. Default:
                                 ${GREEN}$outdir${NOCOLOR}
    -figdir      <figdir>      : Directory of output figures showing model-learned 
                                 brain-behavioral patterns. Default:
                                 ${GREEN}$figdir${NOCOLOR}
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
        -outdir) outdir=$1; shift ;;
        -figdir) figdir=$1; shift ;;
        *) echo "Unknown flag: $flag"; usage; 1>&2; exit 1 
    esac
done

main
