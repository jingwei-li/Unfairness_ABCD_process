#!/bin/sh
# To save space, some files required for re-compute RSFC were deleted. Hence the preprocessing needs to be rerun.

cd ~/storage/ABCD_scripts
source ./config/CBIG_ABCD_paths_t1_fMRI.sh
source ./2020_ABCD_rsfMRI_process/config_rs_GSR_FD0.3_DVARS50/CBIG_ABCD_paths_all.sh
echo $ABCD_table_dir

cd ./2020_ABCD_rsfMRI_process
#./ABCD_rs_proc_wrapper.sh /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt
./ABCD_rs_proc_wrapper.sh /home/jingweil/storage/MyProject/fairAI/ABCD_race/mat/RSFC/no_censor/empty_subjects.txt