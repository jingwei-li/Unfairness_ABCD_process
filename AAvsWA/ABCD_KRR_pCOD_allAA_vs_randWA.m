function ABCD_KRR_pCOD_allAA_vs_randWA(model_dir, bhvr_ls, subj_ls, ...
    split_dir, split_fstem, Nsplits, outmat)

% ABCD_KRR_pCOD_allAA_vs_randWA(model_dir, bhvr_ls, subj_ls, ...
%     split_dir, split_fstem, Nsplits, outmat)
%
% Input:
%   - model_dir
%     The directory storing kernel regression results (full path).
%
%   - bhvr_ls (optional)
%     Behavior list (full path, text file).
%     Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
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
%   - outmat
%     Output filename (full path).

%% default arguments
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, '/behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_ls);

if(~exist('Nsplits', 'var') || isempty(Nsplits))
    Nsplits = nchoosek(10, 3);
end

%% compute predictive COD
pCOD_AA = nan(nbhvr, Nsplits); pCOD_WA = pCOD_AA; pCOD_AAWA = pCOD_AA;
ss_res_AA = pCOD_AA; ss_res_WA = pCOD_WA; ss_res_AAWA = pCOD_AA;
ss_total = pCOD_AA;

AA_pred = cell(nbhvr, Nsplits); WA_pred = AA_pred; AA_test = AA_pred;
WA_test = AA_pred; AA_train = AA_pred; WA_train = AA_pred;

all_selAA = cell(nbhvr,1); all_selWA = all_selAA;
for b = 1:nbhvr
    fprintf('#%d behavior: %s ...\n', b, bhvr_nm{b})
    fold_file = fullfile(split_dir, 'allAA_randWA', [bhvr_nm{b} split_fstem '.mat']);
    load(fold_file)
    if(length(sub_fold) ~= Nsplits)
        error('Nsplits does not equal to length of sub_fold.')
    end
    
    krr_file = fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']);
    krr = load(krr_file);
    
    for f = 1:Nsplits
        all_selAA{b} = [all_selAA{b}; sub_fold(f).selAA];
        all_selWA{b} = [all_selWA{b}; sub_fold(f).selWA];
    end

    for f = 1:Nsplits
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
        [pCOD_AA(b,f), pCOD_WA(b,f), pCOD_AAWA(b,f), ss_res_AA(b,f), ...
            ss_res_WA(b,f), ss_res_AAWA(b,f), ss_total(b,f)] = ...
            ABCD_pCOD_2groups(AA_pred{b,f}, WA_pred{b,f}, AA_test{b,f}, ...
            WA_test{b,f}, AA_train{b,f}, WA_train{b,f});
    end

end

save(outmat, 'pCOD_AA', 'pCOD_WA', 'pCOD_AAWA', 'ss_res_AA', 'ss_res_WA', 'ss_res_AAWA', ...
    'ss_total', 'AA_pred', 'WA_pred', 'AA_test', 'WA_test', 'AA_train', 'WA_train');
    
end