#!/bin/sh
#
# Jingwei Li, 20200820

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"


########################
# input setup
########################
proj_dir=/data/users/jingweil/storage/MyProject/fairAI/ABCD_race

bhvr_ls=$proj_dir/scripts/lists/behavior_list.txt
full_subj_ls=$proj_dir/scripts/lists/subjects_pass_rs_pass_pheno.txt
full_FC_file=$proj_dir/mat/RSFC/pass_rs_pass_pheno_5351.mat
outdir=$proj_dir/models/KRR/20200721/reg_AgeSexMtIcvPEduc_y_randWA
split_dir=$proj_dir/mat/matchANDsplit/20200719/train_randWA
split_fstem="_pass_rs_pass_pheno"

main() {
    mkdir -p $outdir/FC
    behaviors=$(cat $bhvr_ls)
    work_dir=$outdir/logs/HPC
	mkdir -p $work_dir

    for b in $behaviors; do
        if [ ! -f $outdir/FSM/FSM_corr.mat ]; then
		    memory=12
	    else
		    memory=8
	    fi

        subj_ls=$split_dir/subj_${b}${split_fstem}.txt
        FC_file=$outdir/FC/$b.mat
        subfold_f=$split_dir/${b}${split_fstem}.mat

	    mkdir -p $outdir/$b
        cmd="$DIR/ABCD_KRR_in_subgroup.sh -bhvr_name $b -full_subj_ls $full_subj_ls "
        cmd="$cmd -subj_ls $subj_ls -subfold_f $subfold_f -full_FC_file $full_FC_file "
	    cmd="$cmd -FC_file $FC_file -outdir $outdir/$b -csvname $csvname -cfds_ls $cfds_ls"
        echo $cmd

        jname=ABCD_KRR_${b}
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime=25:00:00 -mem \
${memory}G -name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out

        sleep 3s
    done
}

#############################
# Function usage
#############################
usage() { GREEN='\033[0;32m'; NOCOLOR='\033[0m'; echo -e "
NAME:
    ABCD_KRR_reg_AgeSexMtIcvPEduc_y_randWA.sh

DESCRIPTION:
    Wrapper script to perform kernel ridge regression only on randomly selected WA 
    (with the same sample size as all AA).
    It calls ABCD_KRR_in_subgroup.sh to first extract functional connectivity of the subset 
    of subjects (i.e. random WA), then run kernel ridge regression only within these subjects.

ARGUMENTS:
    -h                            : Display help information.
    -csvname       <csvname>      : Full path of the csv file containing all phenotypical values which are 
                                    necessary for this study. 
                                    It is generated by ../preparation/ABCD_read_all_measures.m. Default:
                                    ${GREEN}$csvname${NOCOLOR}
    -bhvr_ls       <bhvr_ls>      : List of behavioral measures involved in this study (full path). Default:
                                    ${GREEN}$bhvr_ls${NOCOLOR}
    -cfds_ls       <cfds_ls>      : List of confounding variables which need to be regressed out (full path).
                                    Default:
                                    ${GREEN}$cfds_ls${NOCOLOR}
    -full_subj_ls  <full_subj_ls> : List of ALL subjects who have passed all quality controls and had all 
                                    required phenotypes (full path). 
                                    It is generated by ../preparation/ABCD_read_all_measures.m. Default:
                                    ${GREEN}$full_subj_ls${NOCOLOR}
    -full_FC_file  <full_FC_file> : Full path of a .mat file containing a 3D matrix which is the resting-state
                                    functional connectivity of all the subjects in <full_subj_ls>. 
                                    It is the output file of ../preparation/ABCD_check_RSFC_NaN.m. Default:
                                    ${GREEN}$full_FC_file${NOCOLOR}
    -split_dir     <split_dir>    : Full path of the directory containing the split folds of the subpopulation
                                    (i.e. random WA).
                                    It is the output directory of ../match_split/ABCD_create_subfold_randWA.m
                                    Default: 
                                    ${GREEN}$split_dir${NOCOLOR}
    -split_fstem   <split_fstem>  : A string shared across all filenames of split folds. It should be the same 
                                    string as passed into ../match_split/ABCD_create_subfold_randWA.m
                                    Default: ${GREEN}$split_fstem${NOCOLOR}
    -outdir        <outdir>       : Full path of output directory. Default:
                                    ${GREEN}$outdir${NOCOLOR}
" 1>&2; exit 1; }

##########################################
# Parse Arguments 
##########################################
while [[ $# -gt 0 ]]; do
    flag=$1; shift;

    case $flag in
        -h) usage; exit ;;
        -csvname) csvname=$1; shift ;;
        -bhvr_ls) bhvr_ls=$1; shift ;;
        -cfds_ls) cfds_ls=$1; shift ;;
        -full_subj_ls) full_subj_ls=$1; shift ;;
        -full_FC_file) full_FC_file=$1; shift ;;
        -split_dir) split_dir=$1; shift ;;
        -split_fstem) split_fstem; shift ;;
        -outdir) outdir=$1; shift ;;
        *) echo "Unknown flag: $flag"; usage; 1>&2; exit 1 
    esac
done

main
