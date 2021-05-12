#!/bin/sh
#
# Wrapper script to perform kernel ridge regression on whole population.
# Confounding variables are regressed from behavioral scores.
#
# Jingwei Li, 20200721

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"


########################
# input setup
########################
proj_dir=/data/users/jingweil/storage/MyProject/fairAI/ABCD_race

ls_dir=$proj_dir/scripts/lists
bhvr_ls=$ls_dir/behavior_list.txt
subj_ls=$ls_dir/subjects_pass_rs_pass_pheno.txt
csvname=$ls_dir/phenotypes_pass_rs.txt
cfds_ls=$ls_dir/confounds_list.txt

FC_file=$proj_dir/mat/RSFC/pass_rs_pass_pheno_5351.mat
outdir=$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_fr_y

subfold_dir=$proj_dir/mat/matchANDsplit/20200719
subfold_stem=_pass_rs_pass_pheno

#########################
# core function
#########################
main() {
	behaviors=$(cat $bhvr_ls)
	work_dir=$outdir/logs/HPC
	mkdir -p $work_dir

	for b in $behaviors; do
		if [ ! -f $outdir/FSM/FSM_corr.mat ]; then
			memory=20
		else
			memory=12
		fi
	
		echo $b
		subfold_f=$subfold_dir/sub_fold${subfold_stem}_${b}.mat
		cmd="$DIR/ABCD_KRR.sh -bhvr_name $b -subj_ls $subj_ls -subfold_f $subfold_f -FC_file \
$FC_file -outdir $outdir -outstem $b -csvname $csvname -cfds_ls $cfds_ls"

		jname=ABCD_KRR_${b}
		$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 25:00:00 -mem ${memory}G \
-name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
		if [ ! -f $outdir/confounds_age_sex_FD_DVARS_ICV_peduc_avg.mat ]; then
			sleep 3m
		else
			sleep 3s
		fi
	done
}

#############################
# Function usage
#############################
usage() { echo "
NAME:
    ABCD_KRR_reg_AgeSexMtIcvPEduc_from_y.sh

DESCRIPTION:
    Wrapper script to perform kernel ridge regression on whole population.
    Confounding variables are regressed from behavioral scores.

ARGUMENTS:
    -h                            : Display help information.
    -csvname       <csvname>      : Full path of the csv file containing all phenotypical values which are 
                                    necessary for this study. 
                                    It is generated by ../preparation/ABCD_read_all_measures.m. Default:
                                    $csvname
    -bhvr_ls       <bhvr_ls>      : List of behavioral measures involved in this study (full path). Default:
                                    $bhvr_ls
    -cfds_ls       <cfds_ls>      : List of confounding variables which need to be regressed out (full path).
                                    Default:
                                    $cfds_ls
    -subj_ls       <subj_ls>      : List of subjects who have passed all quality controls and had all required 
                                    phenotypes (full path). 
                                    It is generated by ../preparation/ABCD_read_all_measures.m. Default:
                                    $subj_ls
    -FC_file       <FC_file>      : Full path of a .mat file containing a 3D matrix which is the resting-state
                                    functional connectivity of all subjects. 
                                    It is the output file of ../preparation/ABCD_check_RSFC_NaN.m. Default:
                                    $FC_file
    -subfold_dir   <subfold_dir>  : The directory containing the split folds. This should be the output 
                                    directory of ../match_split/ABCD_match_and_split.m. Default:
                                    $subfold_dir
    -subfold_stem  <subfold_stem> : A string shared across all filenames of split folds. It should be the same 
                                    string as passed into ../match_split/ABCD_match_and_split.m, 'outstem' 
                                    argument. Default: $subfold_stem
    -outdir        <outdir>       : Full path of output directory. Default:
                                    $outdir

EXAMPLE:
    $DIR/ABCD_KRR_reg_AgeSexMtIcvPEduc_from_y.sh \\
    -csvname <path_to_lists>/phenotypes_pass_rs.txt -bhvr_ls <path_to_lists>/behavior_list.txt -cfds_ls \\
    <path_to_lists>/confounds_list.txt -subj_ls <path_to_lists>/subjects_pass_rs_pass_pheno.txt \\
    -FC_file <path_to_RSFC>/pass_rs_pass_pheno_5351.mat -subfold_dir /.../matchANDsplit/ -subfold_stem \\
    _pass_rs_pass_pheno -outdir /.../reg_AgeSexMtIcvPEduc_fr_y/
" 1>&2; exit 1; }


##########################################
# Parse Arguments 
##########################################

while [[ $# -gt 0 ]]; do
    flag=$1; shift;

    case $flag in
		-h) usage; exit ;;
		-csvname)
			csvname=$1; shift ;;
		-bhvr_ls)
			bhvr_ls=$1; shift ;;
		-cfds_ls)
			cfds_ls=$1; shift ;;
		-subj_ls)
			subj_ls=$1; shift ;;
		-FC_file)
			FC_file=$1; shift ;;
		-subfold_dir)
			subfold_dir=$1; shift ;;
		-subfold_stem)
			subfold_stem=$1; shift ;;
		-outdir)
			outdir=$1; shift ;;
		*)
            echo "Unknown flag: $flag"
            usage; 1>&2; exit 1;;
	esac
done

main
