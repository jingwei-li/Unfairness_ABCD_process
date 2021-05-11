function [FD, DVARS] = ABCD_read_motion(fmri_dir, subj_list, dosave, outdir, outstem)

% [FD, DVARS] = ABCD_read_motion(fmri_dir, subj_list, dosave, outdir, outstem)
%
% Collect the head motion measures of all subjects from the CBIG preprocessed rs-fMRI data.
%
% Inputs:
% - fmri_dir
%   Full path of the directory storing the preprocessed resting-state fMRI
%   data using the CBIG preprocessing pipeline. Measures of head motion 
%   need to be read from this directory.
%   1. FD: [fmri_dir '/' subj_ID '/bold/mc/' subj_ID '_bld' run_ID '_rest_mc_motion_outliers_FDRMS']
%   2. DVARS: [fmri_dir '/' subj_ID '/bold/mc/' subjects{s} '_bld' run_ID '_rest_mc_motion_outliers_DVARS']
%   Default: '/mnt/isilon/CSC2/Yeolab/Data/ABCD/process/y0/rs_GSR'
%
% - subj_list
%   List of subjects which passed fMRI prepreocessing quality control (full path). Default:
%   '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt'
% 
% - dosave
%   A 1/0 value to indicate whether the output text files will be saved or not. 
%
% - outdir
%   Full path of output directory.
% 
% - outstem
%   A shared string attached to all output files. There will be 5 output files:
%   1. [outdir '/FD' outstem '.txt']   -- collection of all subjects' FD
%   2. [outdir '/DV' outstem '.txt']   -- collection of all subjects' DVARS
%   3. [outdir '/subjects' outstem '_pass_pheno.txt']  
%      -- Subject IDs with non-empty records for any of the required measures
%   4. [outdir '/behavior_list.txt']   -- summary of all required behavioral names
%   5. [outdir '/colloquial_list.txt']  -- summary of colloquial names of the behavioral measures
%
% Outputs:
% - FD
%   A #subjects x 1 vector. Each entry is the framewise displacement of one subject.
%
% - DVARS
%   A #subjects x 1 vector. Each entry is the DVARS of one subject.
%
% Author: Jingwei Li


if(~exist('fmri_dir', 'var') || isempty(fmri_dir))
    fmri_dir = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/process/y0/rs_GSR';
end

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt';
end

[subjects, nsub] = CBIG_text2cell(subj_list);

FD = zeros(nsub, 1);
DVARS = zeros(nsub, 1);
for s = 1:nsub
    run_list = fullfile(fmri_dir, subjects{s}, 'logs', [subjects{s} '.bold']);
    
    fid = fopen(run_list);
    runs = fgetl(fid);
    fclose(fid);
    runs = strsplit(runs, ' ');
    empty_idx = cellfun(@isempty, runs);
    runs(empty_idx) = [];
    
    curr_FD = zeros(length(runs), 1);
    curr_DVARS = zeros(length(runs), 1);
    
    for r = 1:length(runs)
        FD_file = fullfile(fmri_dir, subjects{s}, 'bold', 'mc', ...
            [subjects{s} '_bld' runs{r} '_rest_mc_motion_outliers_FDRMS']);
        DV_file = fullfile(fmri_dir, subjects{s}, 'bold', 'mc', ...
            [subjects{s} '_bld' runs{r} '_rest_mc_motion_outliers_DVARS']);
        
        curr_FD(r) = mean(dlmread(FD_file));
        curr_DVARS(r) = mean(dlmread(DV_file));
    end
    FD(s) = mean(curr_FD);
    DVARS(s) = mean(curr_DVARS);
end

if(dosave==1)
    if(~exist(outdir, 'dir'))
        mkdir(outdir)
    end
    FD_save = fullfile(outdir, ['FD' outstem '.txt']);
    DV_save = fullfile(outdir, ['DV' outstem '.txt']);
    
    dlmwrite(FD_save, FD);
    dlmwrite(DV_save, DVARS);
end

end

