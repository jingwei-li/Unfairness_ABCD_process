#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir="/home/jingweil/storage/MyProject/fairAI"
outdir="$proj_dir/cross_ABCD_HCP/models/KRR/no_reg/ABCD_to_HCP"
ABCD_KRR_dir="$proj_dir/ABCD_race/models/KRR/20200721/no_reg"
ABCD_split_dir="$proj_dir/ABCD_race/mat/matchANDsplit/20200719"
crs_dt_FSM_dir="$proj_dir/cross_ABCD_HCP/mat/FSM/no_reg/ABCD_to_HCP"
HCP_KRR_dir="$proj_dir/HCP_race/trained_model/split_948sub_AA_matchedWA_rm_AA_outliers18/outputs/l2_0_20_opt_pCOD_no_reg"
max_HCP_seed=400
HCP_splitAAdir="$proj_dir/HCP_race/mat/split_AA_948_rm_outliers18"
HCP_splitWAdir="$proj_dir/HCP_race/mat/split_WA_rm_AA_outliers18"
usable_seeds_lsdir="$HCP_splitWAdir/usable_seeds"
HCP_subj_ls="$proj_dir/HCP_race/scripts/lists/subjects_wIncome_948.txt"

HCP_comm_bhvr_ls="$proj_dir/cross_ABCD_HCP/scripts/lists/HCP_behavior.txt"
ABCD_comm_bhvr_ls="$proj_dir/cross_ABCD_HCP/scripts/lists/ABCD_behavior.txt"
HCP_comm_bhvr=($(cat $HCP_comm_bhvr_ls))

log_dir=$outdir/logs
work_dir=$log_dir/HPC
mkdir -p $work_dir
i=0
while IFS= read -r ABCD_bhvr_nm; do
    HCP_bhvr_nm=${HCP_comm_bhvr[$i]}
    LF=$log_dir/${ABCD_bhvr_nm}_${HCP_bhvr_nm}.log
    cmd="matlab -nodesktop -nodisplay -nojvm -r \" addpath $DIR; ABCD_KRR_test_HCP( \
        '$outdir', '$ABCD_KRR_dir', '$ABCD_split_dir', '$crs_dt_FSM_dir', \
        '$ABCD_bhvr_nm', '$HCP_bhvr_nm', '$HCP_KRR_dir', $max_HCP_seed, '$HCP_splitAAdir', \
        '$HCP_splitWAdir', '$usable_seeds_lsdir', '$HCP_subj_ls'); exit; \" "

    echo $cmd
    jname=crs_dt_KRR_${ABCD_bhvr_nm}
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 04:00:00 -mem 70G -ncpus 1 \
        -name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
    i=$((i + 1))
done < $ABCD_comm_bhvr_ls