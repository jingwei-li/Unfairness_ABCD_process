function ABCD_fake_discard_frame_lists(outdir, ABCD_rerun_dir, subj_ls)

% ABCD_fake_discard_frame_lists(ourdir, ABCD_rerun_dir, subj_ls)
%
%

if(~exist('ABCD_rerun_dir', 'var') || isempty(ABCD_dir))
    ABCD_rerun_dir = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/process/y0/rs_GSR_test_FD0.3_DVARS50/Jingwei_revision';
end
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt';
end

[subjects, nsub] = CBIG_text2cell(subj_ls);

%% step1: pick the (subject, run) with shortest frames
N_keep = 375;
for s = 1:nsub
    for r = 1:4
        censor_ls = fullfile(ABCD_rerun_dir, subjects{s}, 'qc', ...
            [subjects{s} '_bld00' num2str(r) '_FDRMS0.3_DVARS50_motion_outliers.txt']);
        fMRI_file = fullfile(ABCD_rerun_dir, subjects{s}, 'bold', ['00' num2str(r)], ...
            [subjects{s} '_bld00' num2str(r) '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08.nii.gz']);

        if(exist(fMRI_file))
            outliers = dlmread(censor_ls);
            if(length(outliers) >= 188)
                N_keep_new = length(find(outliers == 1));
                if(N_keep_new < N_keep)
                    N_keep = N_keep_new;
                    shortest_subj = subjects{s};
                    shortest_run = ['00' num2str(r)];
                end
            end
        end
    end
end

fprintf('Shortest subject: %s, run %s, with %d frames\n', shortest_subj, shortest_run, N_keep)

%% step2: truncate every run to be the same length as N_keep
for s = 1:nsub
    for r = 1:4
        censor_ls = fullfile(ABCD_rerun_dir, subjects{s}, 'qc', ...
            [subjects{s} '_bld00' num2str(r) '_FDRMS0.3_DVARS50_motion_outliers.txt']);
        fMRI_file = fullfile(ABCD_rerun_dir, subjects{s}, 'bold', ['00' num2str(r)], ...
            [subjects{s} '_bld00' num2str(r) '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08.nii.gz']);

        if(exist(fMRI_file))
            outliers = dlmread(censor_ls);
            if(length(outliers) >= 188)
                N_keep_new = length(find(outliers == 1));
                outtxt = fullfile(outdir, subjects{s}, ['00' num2str(r) '_truncate_censoring.txt']);
                if(~exist(fullfile(outdir, subjects{s}), 'dir'))
                    mkdir(fullfile(outdir, subjects{s}))
                end
                if(N_keep_new > N_keep)
                    idx = find(outliers==1);
                    len_diff = N_keep_new - N_keep;   % how many frames need to be truncated
                    %fprintf('%s %s, N_keep: %d, N_keep_new: %d, trunctate from: %d\n', ...
                    %    subjects{s}, ['00' num2str(r)], N_keep, N_keep_new, idx(length(idx)-len_diff+1))
                    outliers_new = outliers;
                    outliers_new(idx( (length(idx)-len_diff+1):end )) = 0;
                    dlmwrite(outtxt, outliers_new);
                else
                    system(sprintf('rsync -avz %s %s', censor_ls, outtxt));
                end
            end
        end
    end
end
    
end