function ABCD_LRR_pCOD_AAvsWA(model_dir, bhvr_ls, colloq_ls, subj_ls, split_dir, split_fstem, Nsplits, ...
    predictable_stats, r_th_ls, outmat, bin_flag)

% ABCD_LRR_pCOD_AAvsWA(model_dir, bhvr_ls, colloq_ls, subj_ls, split_dir, split_fstem, Nsplits, ...
%     predictable_stats, outmat)
%
% Calculate predictive COD of matched AA and WA for the linear ridge regression models trained on 
% whole population.
%
% Input:
%   - model_dir
%     The directory storing linear regression results (full path).
%
%   - bhvr_ls (optional)
%     Behavior list (full path, text file).
%     Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
% 
%   - colloq_ls (optional)
%     List of behaviors' colloquial name (full path).
%     Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt'
%
%   - subj_ls (optional)
%     Subject list (full path). Default: 
%  '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt'
%
%   - split_dir
%     The directory storing the data split (full path). It contains a file
%     called ['sel_AAWA' split_fstem '.mat']. For each behaviors, there
%     should be a file called ['sub_fold' split_fstem '_' <behavior name>
%     '.mat'].
%
%   - split_fstem
%     The string that was attached to the filenames of the data split
%     files.
%
%   - Nsplits (optional)
%     Number of splits. Default: 120 (10 choose 3).
%  
%   - predictable_stats
%     Filename of the permutation test result of predictability (full
%     path).
%
%   - r_th_ls
%     Filename of the list of behaviors whose correlation accuracy across all participants is above a threshold.
%
%   - outmat
%     Output filename (full path).
%
%   - bin_flag
%     A 1/0 value, whether the behavioral measure to be predicted is binary or not.
%
% Author: Jingwei Li

%% default arguments
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, '/behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
colloq_nm = CBIG_text2cell(colloq_ls);

if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_ls);

if(~exist('Nsplits', 'var') || isempty(Nsplits))
    Nsplits = nchoosek(10, 3);
end


%% load all seletect AA, WA
match = load(fullfile(split_dir, ['sel_AAWA' split_fstem '.mat']));
all_selAA = cell(nbhvr,1); all_selWA = all_selAA;
for b = 1:nbhvr
    all_selAA{b} = cat(2, match.selAA{b,:});
    all_selWA{b} = cat(1, match.selWA{b,:});
end

%% compute predictive COD
pCOD_AA = nan(nbhvr, Nsplits); pCOD_WA = pCOD_AA; pCOD_AAWA = pCOD_AA;
ss_res_AA = pCOD_AA; ss_res_WA = pCOD_WA; ss_res_AAWA = pCOD_AA;
ss_total = pCOD_AA;

AA_pred = cell(nbhvr, Nsplits); WA_pred = AA_pred; AA_test = AA_pred;
WA_test = AA_pred; AA_train = AA_pred; WA_train = AA_pred;
for b = 1:nbhvr
    fprintf('#%d behavior: %s ...\n', b, bhvr_nm{b})
    fold_file = fullfile(split_dir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']);
    load(fold_file)
    if(length(sub_fold) ~= Nsplits)
        error('Nsplits does not equal to length of sub_fold.')
    end
    
    opt_file = fullfile(model_dir, bhvr_nm{b}, 'results', 'optimal_acc', [bhvr_nm{b} '.mat']);
    opt = load(opt_file);
    
    for f = 1:length(sub_fold)
        lrry = load(fullfile(model_dir, bhvr_nm{b}, 'y', ['fold_' num2str(f)], ...
            ['y_regress_' bhvr_nm{b} '.mat']));
        
        %% collect true & predicted scores of test AA or WA subjects
        AAidx_in_test = zeros(length(sub_fold(f).subject_list), 1);
        WAidx_in_test = zeros(length(sub_fold(f).subject_list), 1);
        [~, idx] = intersect(sub_fold(f).subject_list, sub_fold(f).selAA, 'stable');
        AAidx_in_test(idx) = 1;
        [~, idx] = intersect(sub_fold(f).subject_list, sub_fold(f).selWA, 'stable');
        WAidx_in_test(idx) = 1;

        AAidx_in_all = zeros(length(sub_fold(f).fold_index), 1);
        WAidx_in_all = zeros(length(sub_fold(f).fold_index), 1);
        [~, idx] = intersect(subjects, sub_fold(f).selAA, 'stable');
        AAidx_in_all(idx) = 1;
        [~, idx] = intersect(subjects, sub_fold(f).selWA, 'stable');
        WAidx_in_all(idx) = 1;

        AA_pred{b,f} = opt.y_predict{f}(logical(AAidx_in_test));
        WA_pred{b,f} = opt.y_predict{f}(logical(WAidx_in_test));
        
        AA_test{b,f} = lrry.y_resid(logical(AAidx_in_all));
        WA_test{b,f} = lrry.y_resid(logical(WAidx_in_all));
        
        %% collect true scores of training AA or WA subjects
        trainAA = setdiff(all_selAA{b}, sub_fold(f).selAA);
        trainWA = setdiff(all_selWA{b}, sub_fold(f).selWA);
        
        [~, idx] = intersect(subjects, trainAA, 'stable');
        AA_train{b,f} = lrry.y_resid(idx);
        [~, idx] = intersect(subjects, trainWA, 'stable');
        WA_train{b,f} = lrry.y_resid(idx);
        
        %% compute predictive COD
        [pCOD_AA(b,f), pCOD_WA(b,f), pCOD_AAWA(b,f), ss_res_AA(b,f), ...
            ss_res_WA(b,f), ss_res_AAWA(b,f), ss_total(b,f)] = ...
            ABCD_pCOD_2groups(AA_pred{b,f}, WA_pred{b,f}, AA_test{b,f}, ...
            WA_test{b,f}, AA_train{b,f}, WA_train{b,f});
    end
end

[outdir] = fileparts(outmat);
mkdir(outdir)
save(outmat, 'pCOD_AA', 'pCOD_WA', 'pCOD_AAWA', 'ss_res_AA', 'ss_res_WA', 'ss_res_AAWA', ...
    'ss_total', 'AA_pred', 'WA_pred', 'AA_test', 'WA_test', 'AA_train', 'WA_train');

%% selecte behaviors with significant predictive COD,
% and positive in either AA or WA
% and whose correlation accuracy is above a threshold (e.g. 0.15)
stats = load(predictable_stats);
sig_idx = zeros(nbhvr,1);
sig_idx(stats.sig_perm_idx) = 1;
idx = sig_idx==1 & (mean(pCOD_AA,2)>0 | mean(pCOD_WA,2)>0);
predictable = bhvr_nm(idx);
r_th_bhvr = CBIG_text2cell(r_th_ls);
[~, rth_idx] = intersect(predictable, r_th_bhvr, 'stable');
predictable = predictable(rth_idx);
if(~exist(fullfile(model_dir, 'lists'), 'dir'))
    mkdir(fullfile(model_dir, 'lists'))
end
CBIG_cell2text(predictable, fullfile(model_dir, 'lists', 'pCOD_predictable.txt'))
predictable = colloq_nm(idx);
predictable = predictable(rth_idx);
CBIG_cell2text(predictable, fullfile(model_dir, 'lists', 'pCOD_predictable_colloquial.txt'))

[outdir, outname, outext] = fileparts(outmat);
pCOD_AA = pCOD_AA(idx,:);  pCOD_WA = pCOD_WA(idx,:);  pCOD_AAWA = pCOD_AAWA(idx,:);
ss_res_AA = ss_res_AA(idx,:);  ss_res_WA = ss_res_WA(idx,:);  ss_res_AAWA = ss_res_AAWA(idx,:);
ss_total = ss_total(idx,:);  AA_pred = AA_pred(idx,:);  WA_pred = WA_pred(idx,:);
AA_test = AA_test(idx,:);  WA_test = WA_test(idx,:);
AA_train = AA_train(idx,:);  WA_train = WA_train(idx,:);

pCOD_AA = pCOD_AA(rth_idx,:);  pCOD_WA = pCOD_WA(rth_idx,:);  pCOD_AAWA = pCOD_AAWA(rth_idx,:);
ss_res_AA = ss_res_AA(rth_idx,:);  ss_res_WA = ss_res_WA(rth_idx,:);  ss_res_AAWA = ss_res_AAWA(rth_idx,:);
ss_total = ss_total(rth_idx,:);  AA_pred = AA_pred(rth_idx,:);  WA_pred = WA_pred(rth_idx,:);
AA_test = AA_test(rth_idx,:);  WA_test = WA_test(rth_idx,:);
AA_train = AA_train(rth_idx,:);  WA_train = WA_train(rth_idx,:);
save(fullfile(outdir, [outname '_predictable' outext]), ...
    'pCOD_AA', 'pCOD_WA', 'pCOD_AAWA', 'ss_res_AA', 'ss_res_WA', 'ss_res_AAWA', ...
    'ss_total', 'AA_pred', 'WA_pred', 'AA_test', 'WA_test', 'AA_train', 'WA_train')

end

