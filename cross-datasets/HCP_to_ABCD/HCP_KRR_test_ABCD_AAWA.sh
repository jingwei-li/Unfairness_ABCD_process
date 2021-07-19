#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir="/home/jingweil/storage/MyProject/fairAI"
outdir="$proj_dir/cross_ABCD_HCP/models/KRR/no_reg/HCP_to_ABCD"
HCP_KRR_dir="$proj_dir/HCP_race/trained_model/split_948sub_AA_matchedWA_rm_AA_outliers18/outputs/l2_0_20_opt_pCOD_no_reg"
max_HCP_seed=400
ABCD_selAAWA="$proj_dir/ABCD_race/mat/matchANDsplit/20200719/sel_AAWA_pass_rs_pass_pheno.mat"
HCP_comm_bhvr_ls="$proj_dir/cross_ABCD_HCP/scripts/lists/HCP_behavior.txt"
ABCD_comm_bhvr_ls="$proj_dir/cross_ABCD_HCP/scripts/lists/ABCD_behavior.txt"
crs_dt_FSM_dir="$proj_dir/cross_ABCD_HCP/mat/FSM/no_reg/HCP_to_ABCD"
ABCD_csv="$proj_dir/ABCD_race/scripts/lists/phenotypes_pass_rs.txt"
ABCD_subj_ls="$proj_dir/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt"
ABCD_bhvr_ls="$proj_dir/ABCD_race/scripts/lists/behavior_list.txt"
HCP_splitAAdir="$proj_dir/HCP_race/mat/split_AA_948_rm_outliers18"
HCP_splitWAdir="$proj_dir/HCP_race/mat/split_WA_rm_AA_outliers18"

ABCD_comm_bhvr=($(cat $ABCD_comm_bhvr_ls))
log_dir=$outdir/logs
work_dir=$log_dir/HPC
mkdir -p $work_dir
i=0
while IFS= read -r HCP_bhvr_nm; do
    ABCD_bhvr_nm=${ABCD_comm_bhvr[$i]}
    LF=$log_dir/${HCP_bhvr_nm}_${ABCD_bhvr_nm}.log
    if [ -f $LF ]; then rm $LF; fi
    cmd="matlab -nodesktop -nojvm -nodesktop -r \" addpath $DIR; HCP_KRR_test_ABCD_AAWA('$outdir', '$HCP_KRR_dir', \
        $max_HCP_seed, '$ABCD_selAAWA', '$HCP_bhvr_nm', '$ABCD_bhvr_nm', \
        '$crs_dt_FSM_dir', '$ABCD_csv', '$ABCD_subj_ls', '$ABCD_bhvr_ls', '$HCP_splitAAdir', '$HCP_splitWAdir'); exit; \" >> $LF 2>&1 "

    echo $cmd
    jname=crs_dt_KRR_${HCP_bhvr_nm}
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 04:00:00 -mem 70G -ncpus 1 \
        -name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
    i=$((i + 1))
done < $HCP_comm_bhvr_ls