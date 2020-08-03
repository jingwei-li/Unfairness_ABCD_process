#!/bin/sh
#
# Jingwei Li, 20200721

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

########################
# setup for CIRC cluster
########################
curr_dir=$(pwd)
work_dir=${HOME}/cluster/

echo $curr_dir
echo $work_dir

if [ ! -d $work_dir ]; then
	mkdir -p $work_dir
fi

cd $work_dir


########################
# input setup
########################
proj_dir=/data/users/jingweil/storage/MyProject/fairAI/ABCD_race

ls_dir=$proj_dir/scripts/lists
bhvr_ls=$ls_dir/behavior_list.txt
subj_ls=$ls_dir/subjects_pass_rs_pass_pheno.txt
cfds_ls=NONE

FC_file=$proj_dir/mat/RSFC/pass_rs_pass_pheno_5351_reg_AgeSexMtIcvPEduc.mat
outdir=$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_FC

behaviors=$(cat $bhvr_ls)

for b in $behaviors; do
	if [ -f $outdir/final_result_${b}.mat ]; then
		continue
	fi
	if [ ! -f $outdir/FSM/FSM_corr.mat ]; then
		memory=20
	else
		memory=12
	fi
	
	echo $b
	subfold_f=$proj_dir/mat/matchANDsplit/20200719/sub_fold_pass_rs_pass_pheno_${b}.mat
	cmd="$DIR/ABCD_KRR.sh -bhvr_name $b -subj_ls $subj_ls -cfds_ls $cfds_ls -subfold_f $subfold_f -FC_file $FC_file -outdir $outdir -outstem $b"
	echo $cmd | $CBIG_SCHEDULER_DIR/qsub -V -q circ-spool -l walltime=28:00:00,mem=${memory}GB -m ae -N ABCD_KRR_${b}
	if [ ! -f $outdir/confounds_none.mat ]; then
		sleep 3m
	else
		sleep 3s
	fi
done
