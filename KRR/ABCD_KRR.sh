#!/bin/sh
#
# General bash wrapper to run kernel ridge regression using the ABCD dataset.
# It directly calls the matlab wrapper: ABCD_KRR.m
#
# Jingwei Li, 20200721

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
ls_dir=$HOME/storage/MyProject/fairAI/ABCD_race/scripts/lists
csvname=$ls_dir/phenotypes_pass_rs.txt
cfds_ls=$ls_dir/confounds_list.txt
cfds_X_ls=NONE
LITE=0
N_inner_folds=10
outstem=

main() {
  mkdir -p $outdir/logs
  LF="$outdir/logs/randseed_${outstem}.log"
  if [ -f $LF ]; then rm $LF; fi
  
  echo "bhvr_name = $bhvr_name" >> $LF
  echo "subj_ls = $subj_ls" >> $LF
  echo "subfold_f = $subfold_f" >> $LF
  echo "FC_file = $FC_file" >> $LF
  echo "outdir = $outdir" >> $LF
  echo "csvname = $csvname" >> $LF
  echo "cfds_ls = $cfds_ls" >> $LF
  echo "cfds_X_ls = $cfds_X_ls" >> $LF
  echo "N_inner_folds = $N_inner_folds" >> $LF
  echo "outstem = $outstem" >> $LF
  
  ############ Call matlab function
  matlab -nodesktop -nosplash -nodisplay -r "addpath $DIR; ABCD_KRR('$csvname', '$bhvr_name', '$cfds_ls', '$cfds_X_ls', \
    '$subj_ls', '$subfold_f', '$FC_file', $N_inner_folds, '$outdir', '$outstem', $LITE); exit;" >> $LF 2>&1

  #if [ "$outstem" != "" ]; then stem="_$outstem"; fi
  #if [ -f $outdir/final_result$stem.mat ]; then
  #  rm $outdir/FSM_innerloop/fold_*/FSM*.mat
  #  rm $outdir/FSM_test/fold_*/FSM*.mat
  #fi
}


#############################
# Function usage
#############################
usage() { echo "
NAME:
    ABCD_KRR.sh
	
DESCRIPTION:
    This script calls the matlab wrapper ABCD_KRR.m to perform kernel regression for one behavior.
	
REQUIRED ARGUMENTS:
    -bhvr_name	    <bhvr_name> : Behavior name, the corresponding header in the csv file.
    -subj_ls	    <subj_ls>   : Subject list, full path.
    -subfold_f	    <subfold_f> : Filename of split folds, full path.
    -FC_file	    <FC_file>   : Functional connectivity filename, full path.
    -outdir         <outdir>    : Output directory, full path.

OPTIONAL ARGUMENTS:
    -csvname        <csvname>   : Full path of the csv file containing all phenotypical values which are 
                                  necessary for this study. 
                                  It is generated by ../preparation/ABCD_read_all_measures.m. Default:
                                  $csvname
    -cfds_ls        <cfds_ls>   : List of confounding variables which need to be regressed out from 
                                  behavioral measures (full path). Default:
                                  $cfds_ls
    -cfds_X_ls      <cfds_X_ls> : List of confounding variables which need to be regressed out from RSFC
                                  (full path). Default: 'NONE'.
    -N_inner_folds  <N>         : The number of hyperparameter-selection cross-validation folds. Default: 10.
    -outstem        <outstem>   : A discrimitive string to be attached to output filenames. Default is empty.

EXAMPLE:
    $DIR/ABCD_KRR.sh \\
    -bhvr_name nihtbx_reading_uncorrected -subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt \\
    -subfold_f <path_to_split_folds>/sub_fold_pass_rs_pass_pheno_nihtbx_reading_uncorrected.mat \\
    -FC_file <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -outdir /.../reg_AgeSexMtIcvPEduc_fr_y/ 
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
		-csvname)     # optional
			csvname=$1
			shift;;
		-bhvr_name)
			bhvr_name=$1
			shift;;
		-cfds_ls)     # optional
			cfds_ls=$1
			shift;;
        -cfds_X_ls)   # optional
            cfds_X_ls=$1
            shift;;
		-subj_ls)
			subj_ls=$1
			shift;;
		-subfold_f)
			subfold_f=$1
			shift;;
		-FC_file)
			FC_file=$1
			shift;;
		-N_inner_folds)     # optional
			N_inner_folds=$1
			shift;;
		-outdir)
			outdir=$1
			shift;;
		-outstem)     # optional
			outstem=$1
			shift;;
        -LITE) LITE=$1; shift;;
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


###############################
# check parameters
###############################
if [ -z "$bhvr_name" ]; then
	arg1err "-bhvr_name"
fi
if [ -z "$subj_ls" ]; then
	arg1err "-subj_ls"
fi
if [ -z "$subfold_f" ]; then
	arg1err "-subjfold_f"
fi
if [ -z "$FC_file" ]; then
	arg1err "-FC_file"
fi
if [ -z "$outdir" ]; then
	arg1err "-outdir"
fi


main


