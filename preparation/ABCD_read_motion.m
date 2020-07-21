function [FD, DVARS] = ABCD_read_motion(fmri_dir, subj_list, dosave, out_dir, out_fstem)

if(~exist('fmri_dir', 'var') || isempty(fmri_dir))
    fmri_dir = '/mnt/eql/yeo13/data/ABCD/rs_GSR';
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
    if(~exist(out_dir, 'dir'))
        mkdir(out_dir)
    end
    FD_save = fullfile(out_dir, ['FD' out_fstem '.txt']);
    DV_save = fullfile(out_dir, ['DV' out_fstem '.txt']);
    
    dlmwrite(FD_save, FD);
    dlmwrite(DV_save, DVARS);
end

end

