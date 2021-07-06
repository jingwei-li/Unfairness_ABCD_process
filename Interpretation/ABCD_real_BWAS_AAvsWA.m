function ABCD_real_BWAS_AAvsWA(model_dir, has_subdir, split_dir, split_fstem, outmat, fig_dir, ...
    full_FC, subj_ls, bhvr_ls, colloq_ls, csvname, subj_hdr, nfolds)

% ABCD_real_BWAS_AAvsWA(model_dir, has_subdir, split_dir, split_fstem, outamt, fig_dir, ...
%     full_FC, subj_ls, bhvr_ls, colloq_ls, csvname, subj_hdr, nfolds)
%
% Compute real brain-behavior association in the test fold for matched AA and WA separately.
%
% Inputs:
%   - split_dir
%     Absolut path to the directory storing fold splits.
%   - split_fstem
%     Common filename stem for the fold splits.
%   - outmat
%     Absolute path to the output .mat file.
%   - fig_dir
%     Absolute path to the directory of output figures (If 'None' is given, do not plot the 
%     covariance matrix.
%   - full_FC (optional)
%     Functional connectivity of all subjects in 'subj_ls' (full path). Default:
%     $HOME/storage/MyProject/fairAI/ABCD_race/mat/RSFC/pass_rs_pass_pheno_5351.mat
%   - subj_ls (optional)
%     Subject list (full path). Default: 
%     $HOME/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt
%   - bhvr_ls (optional)
%     Behavior list (full path). Default:
%     $HOME/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt
%   - colloq_ls (optional)
%     List of colloquial names of behaviors (full path). Default:
%     $HOME/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt
%   - csvname (optional)
%     Absolute path to the csv file containing all behavioral measures. Default:
%     $HOME/storage/MyProject/fairAI/ABCD_race/scripts/lists/phenotypes_pass_rs.txt
%   - subj_hdr (optional)
%     Header of the subject ID column in "csvname". Default: 'subjectkey'.
%   - nfolds (optional)
%     Number of fold splits. Default: 120 (10 choose 3).

%% default of optional input arguments
proj_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race');
if(~exist('full_FC', 'var') || isempty(full_FC))
    full_FC = fullfile(proj_dir, 'mat', 'RSFC', 'pass_rs_pass_pheno_5351.mat');
end
load(full_FC)

ls_dir = fullfile(proj_dir, 'scripts', 'lists');
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[full_subj, nsub] = CBIG_text2cell(subj_ls);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
if(exist(bhvr_ls, 'file'))
    [bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
else
    bhvr_nm = {bhvr_ls};  nbhvr = 1;
end

if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
if(exist(colloq_ls, 'file'))
    colloq_nm = CBIG_text2cell(colloq_ls);
else
    colloq_nm = {colloq_ls};
end

if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end
% read csv file
% d = readtable(csvname);
% if(~exist('subj_hdr', 'var') || isempty(subj_hdr))
%     subj_hdr = 'subjectkey';
% end
% allsubj = d.(subj_hdr);
% [~, ~, idx] = intersect(full_subj, allsubj, 'stable');
% d = d(idx, :);

if(~exist('nfolds', 'var') || isempty(nfolds))
    nfolds = nchoosek(10,3);
end

%% colloct y and FC for each fold, each behavior, and compute 
% cov[FC, raw behavior] for test AA and WA separately
cov_testAA = zeros(size(corr_mat, 1), size(corr_mat, 2), nfolds, nbhvr);
cov_testWA = zeros(size(corr_mat, 1), size(corr_mat, 2), nfolds, nbhvr);
for b = 1:nbhvr
    split_fname = fullfile(split_dir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']);
    if(~exist(split_fname, 'file'))
        split_fname = fullfile(split_dir, [bhvr_nm{b} split_fstem '.mat']);
    end
    load(split_fname)
    assert(length(sub_fold) == nfolds, 'nfolds not equal to the length of sub_fold.')

    %y = d.(bhvr_nm{b});
    for f = 1:nfolds
        if(has_subdir)
            yname = fullfile(model_dir, bhvr_nm{b}, 'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']);
        else
            yname = fullfile(model_dir, 'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']);
        end
        y = load(yname);
        y = y.y_resid;

        [~, ~, idxAA] = intersect(sub_fold(f).selAA, full_subj, 'stable');
        y_AA = y(idxAA);
        FC_AA = corr_mat(:,:,idxAA);

        [~, ~, idxWA] = intersect(sub_fold(f).selWA, full_subj, 'stable');
        y_WA = y(idxWA);
        FC_WA = corr_mat(:,:,idxWA);

        cov_testAA(:,:,f,b) = ABCD_cov_FC_behavior(FC_AA, y_AA);
        cov_testWA(:,:,f,b) = ABCD_cov_FC_behavior(FC_WA, y_WA);
    end
end
avg_cov_testAA = squeeze(mean(cov_testAA, 3));
avg_cov_testWA = squeeze(mean(cov_testWA, 3));
%% save
outdir = fileparts(outmat);
mkdir(outdir)
save(outmat, 'cov_testAA', 'avg_cov_testAA', 'cov_testWA', 'avg_cov_testWA', '-v7.3')

%% plot
if(~isempty(fig_dir) && ~strcmpi(fig_dir, 'none'))
    mkdir(fig_dir)
    for b = 1:nbhvr
        CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid(avg_cov_testAA(1:200, 1:200, b), ...
            avg_cov_testAA(1:200, 201:400, b), avg_cov_testAA(201:400, 201:400, b), ...
            avg_cov_testAA(1:200, 401:end, b), avg_cov_testAA(201:400, 401:end, b), ...
            avg_cov_testAA(401:end, 401:end, b), [min(min(min(avg_cov_testAA(:,:,b))), min(min(avg_cov_testWA(:,:,b))))  ...
            max(max(max(avg_cov_testAA(:,:,b))), max(max(avg_cov_testWA(:,:,b))))], fullfile(fig_dir, [colloq_nm{b} '_testAA']))

        CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid(avg_cov_testWA(1:200, 1:200, b), ...
            avg_cov_testWA(1:200, 201:400, b), avg_cov_testWA(201:400, 201:400, b), ...
            avg_cov_testWA(1:200, 401:end, b), avg_cov_testWA(201:400, 401:end, b), ...
            avg_cov_testWA(401:end, 401:end, b), [min(min(min(avg_cov_testAA(:,:,b))), min(min(avg_cov_testWA(:,:,b))))  ...
            max(max(max(avg_cov_testAA(:,:,b))), max(max(avg_cov_testWA(:,:,b))))], fullfile(fig_dir, [colloq_nm{b} '_testWA']))

        CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid(avg_cov_testWA(1:200, 1:200, b) - avg_cov_testAA(1:200, 1:200, b), ...
            avg_cov_testWA(1:200, 201:400, b) - avg_cov_testAA(1:200, 201:400, b), avg_cov_testWA(201:400, 201:400, b) - avg_cov_testAA(201:400, 201:400, b), ...
            avg_cov_testWA(1:200, 401:end, b) - avg_cov_testAA(1:200, 401:end, b), avg_cov_testWA(201:400, 401:end, b) - avg_cov_testAA(201:400, 401:end, b), ...
            avg_cov_testWA(401:end, 401:end, b) - avg_cov_testAA(401:end, 401:end, b), [min(min(avg_cov_testAA(:,:,b) - avg_cov_testWA(:,:,b)))  ...
            max(max(avg_cov_testAA(:,:,b) - avg_cov_testWA(:,:,b)))], fullfile(fig_dir, [colloq_nm{b} '_testAAvsWA']))
    end
end

end