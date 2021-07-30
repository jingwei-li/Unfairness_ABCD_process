#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
Niter=100

proj_dir="/home/jingweil/storage/MyProject/fairAI/ABCD_race"
outdir="$proj_dir/same_size_AAWA/sub_folds"
lsdir="$proj_dir/same_size_AAWA/lists"
orig_split_dir="$proj_dir/mat/matchANDsplit/20200719"

main() {
    work_dir="$outdir/logs"
    mkdir -p $work_dir

    for i in $(seq 2 1 $Niter); do
        cmd="matlab -nodesktop -nojvm -nodisplay -r \"addpath $DIR; ABCD_resample_split('$outdir/sample$i', \
'$lsdir/subjects_sample$i.txt', '$orig_split_dir'); exit; \""
        jname="ABCD_resample_split_$i"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 02:00:00 -mem 2G \
-name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
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
		-orig_split_dir)
			orig_split_dir=$1; shift ;;
		-lsdir)
			lsdir=$1; shift ;;
		-outdir)
			outdir=$1; shift ;;
		*)
            echo "Unknown flag: $flag"
            usage; 1>&2; exit 1;;
	esac
done

main
