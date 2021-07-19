#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir=/home/jingweil/storage/MyProject/fairAI
ABCD_FC=$proj_dir/ABCD_race/mat/RSFC/pass_rs_pass_pheno_5351.mat
HCP_FC=$proj_dir/HCP_race/mat/RSFC_948.mat
ABCD_bhvr_ls=$proj_dir/cross_ABCD_HCP/scripts/lists/ABCD_behavior.txt
HCP_bhvr_ls=$proj_dir/cross_ABCD_HCP/scripts/lists/HCP_behavior.txt
HCP_bhvrs=($(cat $HCP_bhvr_ls))
outdir=$proj_dir/cross_ABCD_HCP/mat/FSM/no_reg/ABCD_to_HCP
work_dir=$outdir/logs
mkdir -p $work_dir

i=0
while IFS= read -r ABCD_bhvr_nm; do
    jname=ABCD_to_HCP_FSM_${ABCD_bhvr_nm}
    HCP_bhvr_nm=${HCP_bhvrs[$i]}
    ABCD_subfold=$proj_dir/ABCD_race/mat/matchANDsplit/20200719/sub_fold_pass_rs_pass_pheno_${ABCD_bhvr_nm}.mat
    cmd="matlab -nodesktop -nojvm -nodisplay -r \" addpath $DIR; ABCD_to_HCP_FSM( \
        '$ABCD_FC', '$HCP_FC', '$ABCD_subfold', '$ABCD_bhvr_nm', '$HCP_bhvr_nm', '$outdir'); exit; \" "
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 08:00:00 -mem 30G -ncpus 1 \
        -name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
    i=$((i + 1))
done < $ABCD_bhvr_ls