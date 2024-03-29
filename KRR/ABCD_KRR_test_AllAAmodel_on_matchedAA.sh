#!/bin/sh
#
# Wrapper bash script to calculate the out-of-sample prediction accuracy of matched AA
# when kernel ridge regression models were trained on all AA.
#
# Author: Jingwei Li

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir="/home/jingweil/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_allAA"

ls_dir=$proj_dir/scripts/lists
full_subj_ls=$ls_dir/subjects_pass_rs_pass_pheno.txt
bhvr_ls="$proj_dir/scripts/lists/behavior_list.txt"
csvname=$ls_dir/phenotypes_pass_rs.txt
cfds_ls=$ls_dir/confounds_list.txt
cfds_X_ls=NONE

full_FC_file=$proj_dir/mat/RSFC/pass_rs_pass_pheno_5351.mat
split_dir="$proj_dir/mat/matchANDsplit/20200719"
split_fstem="_pass_rs_pass_pheno"
AA_subfolder="train_allAA"

main() {
    for b in $(cat $bhvr_ls); do
	    cmd="matlab -nojvm -nodesktop -nodisplay -r \"addpath $DIR; \
ABCD_KRR_test_AllAAmodel_on_matchedAA('$csvname', '$model_dir', '$b', '$cfds_ls', '$cfds_X_ls', \
'$full_subj_ls', '$full_FC_file', '$split_dir', '$split_fstem', '$AA_subfolder'); \
exit;\" "
	    echo $cmd
	    work_dir="$model_dir/$b/logs"
	
	    jname=test_AllAAmodel_on_matchedAA_${b}
	    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 02:00:00 -mem 12G \
-name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
    done
}

#############################
# Function usage
#############################
usage() { GREEN='\033[0;32m'; NOCOLOR='\033[0m'; echo -e "
NAME:
    ABCD_KRR_test_AllAAmodel_on_matchedAA.sh

DESCRIPTION:
    Wrapper bash script to calculate the out-of-sample predicted behavioral scores of matched AA when kernel 
    ridge regression models were trained on all AA.
    The models trained on all AA were generated by ABCD_KRR_reg_AgeSexMtIcvPEduc_y_allAA.sh.

ARGUMENTS:
    -h                            : Display help information.
    -model_dir     <model_dir>    : Full path of the directory containing the kernel ridge regression models 
                                    which were trained on all African Americans. It is the output directory 
                                    of ABCD_KRR_reg_AgeSexMtIcvPEduc_y_allAA.sh. Default:
                                    ${GREEN}$model_dir${NOCOLOR}
    -csvname       <csvname>      : Full path of the csv file containing all phenotypical values which are 
                                    necessary for this study. It is generated by 
                                    ../preparation/ABCD_read_all_measures.m. Default:
                                    ${GREEN}$csvname${NOCOLOR}
    -bhvr_ls       <bhvr_ls>      : List of behavioral measures involved in this study (full path). Default:
                                    ${GREEN}$bhvr_ls${NOCOLOR}
    -cfds_ls       <cfds_ls>      : List of confounding variables which need to be regressed out from 
                                    behavioral scores (full path). Default:
                                    ${GREEN}$cfds_ls${NOCOLOR}
    -cfds_X_ls     <cfds_X_ls>    : List of confounding variables which need to be regressed out from RSFC
                                    (full path). Default: NONE
    -full_subj_ls  <full_subj_ls> : List of ALL subjects who have passed all quality controls and had all 
                                    required phenotypes (full path). It is generated by 
                                    ../preparation/ABCD_read_all_measures.m. Default:
                                    ${GREEN}$full_subj_ls${NOCOLOR}
    -full_FC_file  <full_FC_file> : Full path of a .mat file containing a 3D matrix which is the resting-state
                                    functional connectivity of all the subjects in <full_subj_ls>. It is the 
                                    output file of ../preparation/ABCD_check_RSFC_NaN.m. Default:
                                    ${GREEN}$full_FC_file${NOCOLOR}
    -split_dir     <split_dir>    : Full path of the directory containing the split folds of the whole population.
                                    It is the output directory of ../match_split/ABCD_match_and_split.m. Default: 
                                    ${GREEN}$split_dir${NOCOLOR}
    -split_fstem   <split_fstem>  : A string shared across all filenames of split folds. It should be the same 
                                    string as passed into ../match_split/ABCD_match_and_split.m. Default: 
                                    ${GREEN}$split_fstem${NOCOLOR}
    -AA_subfolder  <AA_subfolder> : The split folds of all AA should have been saved in a subfolder under
                                    <split_dir>. <AA_subfolder> is the relative name of the subfolder. Default: 
                                    ${GREEN}$AA_subfolder${NOCOLOR}
" 1>&2; exit 1; }

##########################################
# Parse Arguments 
##########################################
while [[ $# -gt 0 ]]; do
    flag=$1; shift;

    case $flag in
        -h) usage; exit ;;
        -model_dir) model_dir=$1; shift ;;
        -csvname) csvname=$1; shift ;;
        -bhvr_ls) bhvr_ls=$1; shift ;;
        -cfds_ls) cfds_ls=$1; shift ;;
        -cfds_X_ls) cfds_X_ls=$1; shift;;
        -full_subj_ls) full_subj_ls=$1; shift ;;
        -full_FC_file) full_FC_file=$1; shift ;;
        -split_dir) split_dir=$1; shift ;;
        -split_fstem) split_fstem; shift ;;
        -AA_subfolder) AA_subfolder=$1; shift ;;
        *) echo "Unknown flag: $flag"; usage; 1>&2; exit 1 
    esac
done

main
