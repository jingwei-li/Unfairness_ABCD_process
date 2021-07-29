function ABCD_FC_no_censor_genlists(outdir, ABCD_rerun_dir, subj_ls)

% ABCD_FC_no_censor_genlists(outdir, ABCD_dir, subj_ls)
%
% 

if(~exist('ABCD_dir', 'var') || isempty(ABCD_dir))
    ABCD_dir = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/process/y0/rs_GSR';
end
if(~exist('ABCD_rerun_dir', 'var') || isempty(ABCD_rerun_dir))
    ABCD_rerun_dir = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/process/y0/rs_GSR_test_FD0.3_DVARS50/Jingwei_revision';
end
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt';
end

[subjects, nsub] = CBIG_text2cell(subj_ls);
lh_lines = {}; rh_lines = {}; subcort_lines = {};
empty_subj = {};
for s = 1:nsub
    curr_lh = ''; curr_rh = ''; curr_subcort = '';
    for r = 1:4
        fMRI_file = fullfile(ABCD_rerun_dir, subjects{s}, 'surf', ['lh.' subjects{s} '_bld00' num2str(r) ...
            '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz']);
        if(exist(fMRI_file, 'file'))
            curr_lh = [curr_lh ' ' fMRI_file];
        end
        
        fMRI_file = fullfile(ABCD_rerun_dir, subjects{s}, 'surf', ['rh.' subjects{s} '_bld00' num2str(r) ...
            '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz']);
        if(exist(fMRI_file, 'file'))
            curr_rh = [curr_rh ' ' fMRI_file];
        end

        fMRI_file = fullfile(ABCD_rerun_dir, subjects{s}, 'bold', ['00' num2str(r)], [subjects{s} '_bld00' num2str(r) ...
            '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08.nii.gz']);
        if(exist(fMRI_file, 'file'))
            curr_subcort = [curr_subcort ' ' fMRI_file];
        end
    end
    
    if(~isempty(curr_lh))
        lh_lines = [lh_lines; {curr_lh}];
        if(~exist(fullfile(outdir, 'individuals', subjects{s}), 'dir'))
            mkdir(fullfile(outdir, 'individuals', subjects{s}))
        end
        CBIG_cell2text({curr_lh}, fullfile(outdir, 'individuals', subjects{s}, 'lh_surf_list.txt'))
    end
    if(~isempty(curr_rh))
        rh_lines = [rh_lines; {curr_rh}];
        CBIG_cell2text({curr_rh}, fullfile(outdir, 'individuals', subjects{s}, 'rh_surf_list.txt'))
    end
    if(~isempty(curr_subcort))
        subcort_lines = [subcort_lines; {curr_subcort}];
        CBIG_cell2text({curr_subcort}, fullfile(outdir, 'individuals', subjects{s}, 'volume_list.txt'))
    end

    if(isempty(curr_lh))
        empty_subj = [empty_subj; subjects(s)];
    end
end

CBIG_cell2text(lh_lines, fullfile(outdir, 'lh_surf_list.txt'))
CBIG_cell2text(rh_lines, fullfile(outdir, 'rh_surf_list.txt'))
CBIG_cell2text(subcort_lines, fullfile(outdir, 'volume_list.txt'))
CBIG_cell2text(empty_subj, fullfile(outdir, 'empty_subjects.txt'))
    
end