## Read (behavioral, demographic, morphological) phenotypes from ABCD csv files. Censor out the subjects with none of the required measures. Combine these measures into one csv file.

The user only needs to run the top-level script `ABCD_read_all_measures.m`. This script call other scriot such as `ABCD_read_gender.m`, `ABCD_read_LittleMan.m`, etc. Usage example (matlab):

```matlab
ABCD_read_all_measures('/preprocessed/rs-fMRI/dir/', '/sMRI/recon_all/dir/', ...
    '/list/of/subjects/passed/rs-fMRI/QC.txt', ...
    '/output/dir/', '<a shared string attached on all output files>')
```

##