function ABCD_FC_truncate_genlists(outdir, fake_censor_dir, ABCD_rerun_dir, subj_ls)

% ABCD_FC_truncate_genlists(outdir, fake_censor_dir, ABCD_rerun_dir, subj_ls)
%
%

if(~exist('ABCD_rerun_dir', 'var') || isempty(ABCD_rerun_dir))
    ABCD_rerun_dir = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/process/y0/rs_GSR_test_FD0.3_DVARS50/Jingwei_revision';
end
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt';
end

[subjects, nsub] = CBIG_text2cell(subj_ls);
for s = 1:nsub
    curr_lh = ''; curr_rh = ''; curr_subcort = ''; curr_fake_censor = '';
    for r = 1:4
        fake_censor_txt = fullfile(fake_censor_dir, subjects{s}, ['00' num2str(r) '_truncate_censoring.txt']);
        if(exist(fake_censor_txt, 'file'))
            curr_fake_censor = [curr_fake_censor ' ' fake_censor_txt];

            fMRI_file = fullfile(ABCD_rerun_dir, subjects{s}, 'surf', ['lh.' subjects{s} '_bld00' num2str(r) ...
                '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz']);
            curr_lh = [curr_lh ' ' fMRI_file];
        
            fMRI_file = fullfile(ABCD_rerun_dir, subjects{s}, 'surf', ['rh.' subjects{s} '_bld00' num2str(r) ...
                '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz']);
            curr_rh = [curr_rh ' ' fMRI_file];

            fMRI_file = fullfile(ABCD_rerun_dir, subjects{s}, 'bold', ['00' num2str(r)], [subjects{s} '_bld00' num2str(r) ...
                '_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08.nii.gz']);
            curr_subcort = [curr_subcort ' ' fMRI_file];
        end
        
    end
    
    if(~isempty(curr_lh))
        if(~exist(fullfile(outdir, 'individuals', subjects{s}), 'dir'))
            mkdir(fullfile(outdir, 'individuals', subjects{s}))
        end

        CBIG_cell2text({curr_lh}, fullfile(outdir, 'individuals', subjects{s}, 'lh_surf_list.txt'))
        CBIG_cell2text({curr_rh}, fullfile(outdir, 'individuals', subjects{s}, 'rh_surf_list.txt'))
        CBIG_cell2text({curr_subcort}, fullfile(outdir, 'individuals', subjects{s}, 'volume_list.txt'))
        CBIG_cell2text({curr_fake_censor}, fullfile(outdir, 'individuals', subjects{s}, 'fake_censor_list.txt'))
    end
end
    
end