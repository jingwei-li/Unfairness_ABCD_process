#!/bin/sh

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

ABCD_dir="/mnt/isilon/CSC2/Yeolab/Data/ABCD/process/y0/rs_GSR"
ROI_dir="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label"
lh_ROI="$ROI_dir/lh.Schaefer2018_400Parcels_17Networks_order.annot"
rh_ROI="$ROI_dir/rh.Schaefer2018_400Parcels_17Networks_order.annot"

proj_dir="/home/jingweil/storage/MyProject/fairAI/ABCD_race"
outdir_parent="$proj_dir/mat/RSFC/no_censor"
output_prefix="RSFC_5351_no_censor"

subj_ls="$proj_dir/scripts/lists/subjects_pass_rs_pass_pheno.txt"

main() {
    for s in $(cat $subj_ls); do
        work_dir="$outdir_parent/logs/individuals/$s"
        mkdir -p $work_dir
        LF="$work_dir/FC_no_censor.log"

        outdir="$outdir_parent/individuals/$s"
        lh_surf_ls="$outdir/lh_surf_list.txt"
        rh_surf_ls="$outdir/rh_surf_list.txt"
        subcort_ls="$outdir/volume_list.txt"

        subcort_ROI="$ABCD_dir/$s/FC_metrics/ROIs/${s}.subcortex.19aseg.func.nii.gz"
        
        cmd="matlab -nodesktop -nojvm -nodisplay -r \"addpath $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/utilities;CBIG_preproc_FCmetrics('$lh_ROI','$rh_ROI','$subcort_ROI','$lh_surf_ls','$rh_surf_ls','$subcort_ls','NONE','Pearson_r','$outdir','$output_prefix');exit;\""

        if [ -f $outdir/${output_prefix}_all2all.mat ]; then continue; fi
        jname="ABCD_FC_no_censor"
        echo $cmd
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 03:00:00 -mem 10G \
            -name $jname -joberr $work_dir/$jname.err -jobout $work_dir/$jname.out
        sleep 3s
    done
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

while [[ $# -gt 0 ]]; do
	flag=$1; shift;
	
	case $flag in
        -h)  usage; 1>&2; exit 1;;
        -ABCD_dir) ABCD_dir=$1; shift;;
        -lh_ROI) lh_ROI=$1; shift;;
        -rh_ROI) rh_ROI=$1; shift;;
        -outdir_parent) outdir_parent=$1; shift;;
        -output_prefix) output_prefix=$1; shift;;
        -subj_ls) subj_ls=$1; shift;;
        *)
            echo "Unknown flag: $flag"
            usage; 1>&2; exit 1
			;;
    esac
done


main
