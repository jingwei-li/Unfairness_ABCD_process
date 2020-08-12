#!/bin/sh
#
# Jingwei Li, 20200811

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

########################################
# main commands to be submitted to jobs
########################################
main(){
    mkdir -p $outdir/logs/
    LF="$outdir/logs/create_FC_${bhvr_name}.log"
    if [ ! -f $FC_file ]; then
        matlab -nodesktop -nojvm -nodisplay -r "addpath(fileparts('$DIR')); \
            ABCD_FC_subgroup('$full_subj_ls', '$subj_ls', '$full_FC_file', '$FC_file'); \
            exit;" >> $LF 2>&1
    fi

    $DIR/ABCD_KRR.sh -bhvr_name $bhvr_name -subj_ls $subj_ls -subfold_f $subfold_f -FC_file $FC_file -outdir $outdir -outstem $bhvr_name

    if [ -f $outdir/final_result_${bhvr_name}.mat ]; then
        echo "KRR run successfully. FC file deleted. :)" >> $LF
        rm $FC_file
    fi
}

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
		-bhvr_name)
			bhvr_name=$1
			shift;;
        -full_subj_ls)
            full_subj_ls=$1
            shift;;
		-subj_ls)
			subj_ls=$1
			shift;;
		-subfold_f)
			subfold_f=$1
			shift;;
        -full_FC_file)
            full_FC_file=$1
            shift;;
		-FC_file)
			FC_file=$1
			shift;;
		-outdir)
			outdir=$1
			shift;;
		*) 
			echo "Unknown flag $flag"
			usage; 1>&2; exit 1
            ;;
	esac
done

main