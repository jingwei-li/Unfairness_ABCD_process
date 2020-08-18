#!/bin/sh
#
# Jingwei Li, 20200811

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

bhvr_ls=$proj_dir/scripts/lists/behavior_list.txt
full_subj_ls=$proj_dir/scripts/lists/subjects_pass_rs_pass_pheno.txt
full_FC_file=$proj_dir/mat/RSFC/pass_rs_pass_pheno_5351.mat
outdir=$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_allAA
mkdir -p $outdir/FC

behaviors=$(cat $bhvr_ls)

for b in $behaviors; do
    if [ ! -f $outdir/FSM/FSM_corr.mat ]; then
		memory=12
	else
		memory=8
	fi

    subj_ls=$proj_dir/mat/matchANDsplit/20200719/allAA/subj_${b}_pass_rs_pass_pheno.txt
    FC_file=$outdir/FC/$b.mat
    subfold_f=$proj_dir/mat/matchANDsplit/20200719/allAA/${b}_pass_rs_pass_pheno.mat

	mkdir -p $outdir/$b
    cmd="$DIR/ABCD_KRR_in_subgroup.sh -bhvr_name $b -full_subj_ls $full_subj_ls "
    cmd="$cmd -subj_ls $subj_ls -subfold_f $subfold_f -full_FC_file $full_FC_file "
	cmd="$cmd -FC_file $FC_file -outdir $outdir/$b"
    echo $cmd
    echo $cmd | $CBIG_SCHEDULER_DIR/qsub -V -q circ-spool -l \
        walltime=25:00:00,mem=${memory}GB -m ae -N ABCD_KRR_${b}

    sleep 3s
done