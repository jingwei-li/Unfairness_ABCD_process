#!/bin/sh
#
# The second step to perform multi-level block permutation test for the predictability
# of linear ridge regression (LRR) across all ABCD subjects.
# In this step, LRR is repeated using permuted behavioral scores. The permutation 
# indicied need to be pre-generated using `ABCD_LRR_predictable_behavior_step1.m`.
#
# Author: Jingwei Li

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir=/home/jingweil/storage/MyProject/fairAI/ABCD_race
model_dir=$proj_dir/models/LRR/20210629/reg_AgeSexMtIcvPEduc_fr_y_FC
cov_stem="_age_gender_FD_DVARS_ICV_peduc_avg"
FC_file=$proj_dir/mat/RSFC/pass_rs_pass_pheno_5351.mat
folds_dir=$proj_dir/mat/matchANDsplit/20200719
folds_fstem="_pass_rs_pass_pheno"
bhvr_ls=$proj_dir/scripts/lists/behavior_list.txt
colloq_ls=$proj_dir/scripts/lists/colloquial_list.txt
Nperm=1000

main() {
    bhvr_nm=($(cat $bhvr_ls))

    i=0
    while IFS= read -r colloq_nm; do
        curr_bhvr=${bhvr_nm[$i]}
        log_dir=$model_dir/$curr_bhvr/logs

        for curr_fold in $(seq 1 1 120); do
            curr_out=$model_dir/$curr_bhvr/perm/fold_${curr_fold}.mat
            if [ -f $curr_out ]; then continue; fi
            LF=$log_dir/predictability_step2_${curr_fold}.log
            if [ -f $LF ]; then rm $LF; fi

            echo "model_dir = $model_dir" >> $LF
            echo "FC_file = $FC_file" >> $LF
            echo "folds_dir = $folds_dir" >> $LF
            echo "folds_fstem = $folds_fstem" >> $LF
            echo "bhvr_nm = $curr_bhvr" >> $LF
            echo "colloq_nm = $colloq_nm" >> $LF

            cmd="matlab -nodesktop -nojvm -nodisplay -r \" addpath $DIR; ABCD_LRR_predictable_behavior_step2(\
                '$model_dir', '$cov_stem', '$FC_file', '$folds_dir', '$folds_fstem', $curr_fold, $Nperm, '$curr_bhvr', \
                '$colloq_nm'); exit; \" >> $LF 2>&1 "
            work_dir=$log_dir/HPC
            jname=ABCD_LRR_predictability_step2_${curr_bhvr}_${curr_fold}
            $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 04:00:00 -mem 30G -ncpus 1 \
                -name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
            sleep 10s
        done
        i=$((i + 1))
    done < $colloq_ls
}


#############################
# Function usage
#############################
usage() { echo "
NAME:

" 1>&2; exit 1; }

##########################################
# Parse Arguments 
##########################################
# Display help message if no argument is supplied
if [ $# -eq 0 ]; then
	usage; 1>&2; exit 1
fi

while [[ $# -gt 0 ]]; do
	flag=$1; shift;
	
	case $flag in
        -model_dir)      # optional
            model_dir=$1; shift;;
        -FC_file)        # optional
            FC_file=$1; shift;;
        -folds_dir)      # optional
            folds_dir=$1; shift;;
        -folds_fstem)    # optional
            folds_fstem=$1; shift;;
        -Nperm)          # optional
            Nperm=$1; shift;;
        -bhvr_ls)        # optional
            bhvr_ls=$1; shift;;
        -colloq_ls)      # optional
            colloq_ls=$1; shift;;
		*) 
			echo "Unknown flag $flag"
			usage; 1>&2; exit 1
			;;
	esac
done

##########################################
# ERROR message
##########################################	
arg1err() {
	echo "ERROR: flag $1 requires one argument"
	exit 1
}



main