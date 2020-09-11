#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"


curr_dir=$(pwd)

proj_dir="/data/users/jingweil/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_allAA"
bhvr_ls="$proj_dir/scripts/lists/behavior_list.txt"
split_dir="$proj_dir/mat/matchANDsplit/20200719"
split_fstem="_pass_rs_pass_pheno"
AA_subfolder="train_allAA"

for b in $(cat $bhvr_ls); do
	cmd="matlab -nojvm -nodesktop -nodisplay -r \"addpath $DIR; \
ABCD_KRR_test_AllAAmodel_on_matchedAA([], '$model_dir', \
'$b', [], [], [], '$split_dir', '$split_fstem', '$AA_subfolder'); \
exit;\" "
	echo $cmd
	work_dir="$model_dir/$b/logs"
	cd $work_dir
	echo $work_dir
	echo $cmd | $CBIG_SCHEDULER_DIR/qsub -V -q circ-spool -l \
walltime=02:00:00,mem=12GB -m ae -N test_AllAAmodel_on_matchedAA_${b} 
	sleep 3s
done