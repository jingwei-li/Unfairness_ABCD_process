function ABCD_KRR_corr_allAA_vs_randWA(model_dir, bhvr_ls, subj_ls, split_dir, split_fstem, ...
    Nsplits, outmat)

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

%% compute Pearson's correlation
corr_AA = nan(nbhvr, Nsplits); corr_WA = corr_AA; corr_AAWA = corr_AA;
AA_pred = cell(nbhvr, Nsplits); WA_pred = AA_pred;
AA_test = AA_pred; WA_test = AA_pred;
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

        %% compute Pearson's r
        corr_AA(b,f) = CBIG_corr(AA_pred{b,f}, AA_test{b,f});
        corr_WA(b,f) = CBIG_corr(WA_pred{b,f}, WA_test{b,f});
        corr_AAWA(b,f) = CBIG_corr([AA_pred{b,f}; WA_pred{b,f}], ...
            [AA_test{b,f}; WA_test{b,f}]);
    end
end

save(outmat, 'corr_AA', 'corr_WA', 'corr_AAWA', 'AA_pred', 'WA_pred', ...
    'AA_test', 'WA_test')

end