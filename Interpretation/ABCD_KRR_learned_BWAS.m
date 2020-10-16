function ABCD_KRR_learned_BWAS(model_dir, has_subdir, split_dir, split_fstem, outmat, fig_dir, ...
    full_FC, subj_ls, bhvr_ls, colloq_ls, nfolds)

% ABCD_KRR_learned_BWAS(model_dir, has_subdir, split_dir, split_fstem, outdir, outstem,...
%     subj_ls, bhvr_ls, colloq_ls)
%
% Compute the brain-wide behavior association learned by kernal regression model.
% In detial, calculate the covariance between FC and predicted behaviors in the training
% set.
% 
% Inputs:
%   - model_dir
%     Absolute path to the directory of KRR models.
%   - has_subdir
%     If KRR results are save for each behavior separately in a subfolder, pass 1, otherwise pass 0.
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

%% Default of optional input arguments
if(~exist('nfolds', 'var') || isempty(nfolds))
    nfolds = nchoosek(10,3);
end

proj_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race');
ls_dir = fullfile(proj_dir, 'scripts', 'lists');
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[full_subj, nsub] = CBIG_text2cell(subj_ls);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
colloq_nm = CBIG_text2cell(colloq_ls);

if(~exist('full_FC', 'var') || isempty(full_FC))
    full_FC = fullfile(proj_dir, 'mat', 'RSFC', 'pass_rs_pass_pheno_5351.mat');
end
load(full_FC)

%% collect FC and behavior for the training set of each fold split, and compute covariance
learned_cov = zeros(size(corr_mat,1), size(corr_mat,2), nfolds, nbhvr);
for b = 1:nbhvr
    fprintf('#%d behavior: %s\n', b, bhvr_nm{b})
    %% load fold splits
    split_fname = fullfile(split_dir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']);
    load(split_fname)
    assert(length(sub_fold) == nfolds, 'Input nfolds deviartes from #folds in sub_fold.')

    for f = 1:nfolds
        %% load KRR prediction of training subjects
        if(has_subdir)
            training = fullfile(model_dir, bhvr_nm{b}, 'test_cv', ['fold_' num2str(f)], ...
                ['opt_training_set_' bhvr_nm{b}, '.mat']);
        else
            training = fullfile(model_dir, 'test_cv', ['fold_' num2str(f)], ...
                ['opt_training_set_' bhvr_nm{b}, '.mat']);
        end
        training = load(training);

        %% select FC of training subjects
        FC = corr_mat(:,:, sub_fold(f).fold_index==0);
        curr_y = training.y_p{1};

        %% compute covariance between FC and predicted scores across training subjects
        learned_cov(:,:,f,b) = ABCD_cov_FC_behavior(FC, curr_y);
    end
end
avg_learned_cov = squeeze(mean(learned_cov, 3));
%% save
outdir = fileparts(outmat);
mkdir(outdir)
save(outmat, 'learned_cov', 'avg_learned_cov', '-v7.3')

%% Plot
if(~isempty(fig_dir) && ~strcmpi(fig_dir, 'none'))
    mkdir(fig_dir)
    for b = 1:nbhvr
        CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid(avg_learned_cov(1:200, 1:200, b), ...
            avg_learned_cov(1:200, 201:400, b), avg_learned_cov(201:400, 201:400, b), ...
            avg_learned_cov(1:200, 401:end, b), avg_learned_cov(201:400, 401:end, b), ...
            avg_learned_cov(401:end, 401:end, b), [min(min(avg_learned_cov(:,:,b)))  ...
            max(max(avg_learned_cov(:,:,b)))], fullfile(fig_dir, colloq_nm{b}))
    end
end
    
end