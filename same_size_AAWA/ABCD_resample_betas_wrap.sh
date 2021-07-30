#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

proj_dir="/home/jingweil/storage/MyProject/fairAI/ABCD_race"
outdir="$proj_dir/same_size_AAWA/betas"
bhvr_ls="$proj_dir/scripts/lists/behavior_list.txt"
csvname="$proj_dir/scripts/lists/phenotypes_pass_rs.txt"
cfds_ls="$proj_dir/scripts/lists/confounds_list.txt"
cfds_X_ls=$cfds_ls
orig_subj_ls="$proj_dir/scripts/lists/subjects_pass_rs_pass_pheno.txt"
RSFC_file="$proj_dir/mat/RSFC/pass_rs_pass_pheno_5351.mat"

new_lsdir="$proj_dir/same_size_AAWA/lists"
sub_fold_dir="$proj_dir/same_size_AAWA/sub_folds"
Niter=100

main() {
    behaviors=$(cat $bhvr_ls)

    for i in $(seq 1 1 $Niter); do
        work_dir="$outdir/logs/sample$i"
        mkdir -p $work_dir

        for bhvr_nm in $behaviors; do
            cmd="matlab -nodesktop -nodisplay -nojvm -r \" addpath $DIR; ABCD_resample_betas( \
'$outdir/sample$i/$bhvr_nm', '$bhvr_nm', '$new_lsdir/subjects_sample$i.txt', '$sub_fold_dir/sample$i/$bhvr_nm.mat', \
'$csvname', '$cfds_ls', '$cfds_X_ls', '$orig_subj_ls', '$RSFC_file'); exit; \""

            jname=ABCD_resample_betas_$bhvr_nm
            $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 03:00:00 -mem 35G \
-name $jname -joberr $work_dir/$bhvr_nm.err -jobout $work_dir/$bhvr_nm.out
        done
        sleep 100s
    done
}

#############################
# Function usage
#############################
usage() { echo "
NAME:
    ABCD_resample_split_wrap.sh
" 1>&2; exit 1; }

##########################################
# Parse Arguments 
##########################################

while [[ $# -gt 0 ]]; do
    flag=$1; shift;

    case $flag in
		-h) usage; exit ;;
		-outdir) outdir=$1; shift;;
        -bhvr_ls) bhvr_ls=$1; shift;;
        -csvname) csvname=$1; shift;;
        -cfds_ls) cfds_ls=$1; shift;;
        -cfds_X_ls) cfds_X_ls=$1; shift;;
        -orig_subj_ls) orig_subj_ls=$1; shift;;
        -RSFC_file) RSFC_file=$1; shift;;
        -new_lsdir) new_lsdir=$1; shift;;
        -sub_fold_dir) sub_fold_dir=$1; shift;;
        -Niter) Niter=$1; shift;;
		*)
            echo "Unknown flag: $flag"
            usage; 1>&2; exit 1;;
	esac
done

main
