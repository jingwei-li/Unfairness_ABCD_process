function ABCD_KRR_pCOD_AAvsWA(model_dir, bhvr_ls, colloq_ls, subj_ls, split_dir, split_fstem, Nsplits, ...
    predictable_stats, outmat)

% ABCD_KRR_pCOD_AAvsWA(model_dir, bhvr_ls, colloq_ls, subj_ls, split_dir, split_fstem, Nsplits, ...
%     predictable_stats, outmat)
%
% Input:
%   - model_dir
%     The directory storing kernel regression results (full path).
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
%   - outmat
%     Output filename (full path).

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
    
    krr_file = fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']);
    krr = load(krr_file);
    
    for f = 1:length(sub_fold)
        krry = load(fullfile(model_dir, 'y', ['fold_' num2str(f)], ...
            ['y_regress_' bhvr_nm{b} '.mat']));
        
        %% collect true & predicted scores of test AA or WA subjects
        AAidx = zeros(nsub, 1);
        WAidx = zeros(nsub, 1);
        [~, idx] = intersect(subjects, sub_fold(f).selAA, 'stable');
        AAidx(idx) = 1;
        [~, idx] = intersect(subjects, sub_fold(f).selWA, 'stable');
        WAidx(idx) = 1;
        
        AA_pred{b,f} = krr.y_predict_concat(logical(AAidx));
        WA_pred{b,f} = krr.y_predict_concat(logical(WAidx));
        
        AA_test{b,f} = krry.y_resid(logical(AAidx));
        WA_test{b,f} = krry.y_resid(logical(WAidx));
        
        %% collect true scores of training AA or WA subjects
        trainAA = setdiff(all_selAA{b}, sub_fold(f).selAA);
        trainWA = setdiff(all_selWA{b}, sub_fold(f).selWA);
        
        [~, idx] = intersect(subjects, trainAA, 'stable');
        AA_train{b,f} = krry.y_resid(idx);
        [~, idx] = intersect(subjects, trainWA, 'stable');
        WA_train{b,f} = krry.y_resid(idx);
        
        %% compute predictive COD
        ss_res_AA(b,f) = sum((AA_pred{b,f} - AA_test{b,f}).^2, 1) ./ length(AA_test{b,f});
        ss_res_WA(b,f) = sum((WA_pred{b,f} - WA_test{b,f}).^2, 1) ./ length(WA_test{b,f});
        ss_res_AAWA(b,f) = sum( ([AA_pred{b,f};WA_pred{b,f}] - [AA_test{b,f};WA_test{b,f}]).^2, 1 ) ...
            ./ length([AA_test{b,f};WA_test{b,f}]);
        ss_total(b,f) = sum( ([AA_train{b,f}; WA_train{b,f}] - mean([AA_train{b,f};WA_train{b,f}],1)).^2, 1 ) ...
            ./ length([AA_train{b,f};WA_train{b,f}]);
        
        pCOD_AA(b,f) = bsxfun(@minus, 1, ss_res_AA(b,f)./ss_total(b,f));
        pCOD_WA(b,f) = bsxfun(@minus, 1, ss_res_WA(b,f)./ss_total(b,f));
        pCOD_AAWA(b,f) = bsxfun(@minus, 1, ss_res_AAWA(b,f)./ss_total(b,f));
    end
end

save(outmat, 'pCOD_AA', 'pCOD_WA', 'pCOD_AAWA', 'ss_res_AA', 'ss_res_WA', 'ss_res_AAWA', ...
    'ss_total', 'AA_pred', 'WA_pred', 'AA_test', 'WA_test', 'AA_train', 'WA_train');

%% selecte behaviors with significant predictive COD,
% and positive in either AA or WA
stats = load(predictable_stats);
sig_idx = zeros(nbhvr,1);
sig_idx(stats.sig_perm_idx) = 1;
idx = sig_idx==1 & (mean(pCOD_AA,2)>0 | mean(pCOD_WA,2)>0);
predictable = bhvr_nm(idx);
if(~exist(fullfile(model_dir, 'lists'), 'dir'))
    mkdir(fullfile(model_dir, 'lists'))
end
CBIG_cell2text(predictable, fullfile(model_dir, 'lists', 'pCOD_predictable.txt'))
predictable = colloq_nm(idx);
CBIG_cell2text(predictable, fullfile(model_dir, 'lists', 'pCOD_predictable_colloquial.txt'))

[outdir, outname, outext] = fileparts(outmat);
pCOD_AA = pCOD_AA(idx,:);  pCOD_WA = pCOD_WA(idx,:);  pCOD_AAWA = pCOD_AAWA(idx,:);
ss_res_AA = ss_res_AA(idx,:);  ss_res_WA = ss_res_WA(idx,:);  ss_res_AAWA = ss_res_AAWA(idx,:);
ss_total = ss_total(idx,:);  AA_pred = AA_pred(idx,:);  WA_pred = WA_pred(idx,:);
AA_test = AA_test(idx,:);  WA_test = WA_test(idx,:);
AA_train = AA_train(idx,:);  WA_train = WA_train(idx,:);
save(fullfile(outdir, [outname '_predictable' outext]), ...
    'pCOD_AA', 'pCOD_WA', 'pCOD_AAWA', 'ss_res_AA', 'ss_res_WA', 'ss_res_AAWA', ...
    'ss_total', 'AA_pred', 'WA_pred', 'AA_test', 'WA_test', 'AA_train', 'WA_train')

end

