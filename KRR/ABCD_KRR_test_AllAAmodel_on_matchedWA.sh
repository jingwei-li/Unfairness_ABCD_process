#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

curr_dir=$(pwd)

proj_dir="/home/jingweil/storage/MyProject/fairAI/ABCD_race"
model_dir="$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_allAA"
bhvr_ls="$proj_dir/scripts/lists/behavior_list.txt"
split_dir="$proj_dir/mat/matchANDsplit/20200719"
split_fstem="_pass_rs_pass_pheno"
AA_subdir="train_allAA"
AAWA_subdir=""
outstem="_matchedWA"

for b in $(cat $bhvr_ls); do
	cmd="matlab -nojvm -nodesktop -nodisplay -r \"addpath $DIR; \
ABCD_KRR_test_AAmodel_on_WA([], '$model_dir', '$b', [], [], [], \
'$split_dir', '$split_fstem', '$AA_subdir', '$AAWA_subdir', '$outstem'); \
exit;\" "
	echo $cmd
	work_dir="$model_dir/$b/logs"
	#cd $work_dir
	jname=test_AllAAmodel_on_matchedWA_${b} 
	$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 02:00:00 -mem 12G \
-name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
	#echo $cmd | $CBIG_SCHEDULER_DIR/qsub -V -q circ-spool -l \
#walltime=02:00:00,mem=12GB -m ae -N test_AllAAmodel_on_matchedWA_${b} 
	sleep 3s
done