#!/bin/sh
#
# Jingwei Li, 20200826

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir="/data/users/jingweil/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_allAA_randWA"
bhvr_ls="$proj_dir/scripts/lists/behavior_list.txt"
split_dir="$proj_dir/mat/matchANDsplit/20200719"
split_fstem="_pass_rs_pass_pheno"
train_folds_subdir="train_allAA_randWA"
matched_subdir=""
outstem="_matchedAA"

for b in $(cat $bhvr_ls); do
	cmd="matlab -nojvm -nodesktop -nodisplay -r \"addpath $DIR; \
ABCD_KRR_test_WAmodel_on_AA([], '$model_dir', '$b', [], [], [], \
'$split_dir', '$split_fstem', '$train_folds_subdir', '$AAWA_subdir', '$outstem'); \
exit;\" "
	echo $cmd
	work_dir="$model_dir/$b/logs"
	cd $work_dir
	echo $cmd | $CBIG_SCHEDULER_DIR/qsub -V -q circ-spool -l \
walltime=02:00:00,mem=12GB -m ae -N test_allAArandWAmodel_on_matchedAA_${b} 
	sleep 3s
done