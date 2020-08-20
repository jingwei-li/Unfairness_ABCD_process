function ABCD_create_subfold_randWA(subj_ls, bhvr_ls, split_dir, split_fstem)

% ABCD_create_subfold_randWA(subj_ls, bhvr_ls, split_dir, split_fstem)
%
% This function is based on the results from `ABCD_select_allAA_randWA.m`.
% `ABCD_select_allAA_randWA.m` has selected all AA and random WA with the 
% same sample size (if #WA < #AA for some site, only part of AA subjects 
% are selected to match the size).
% The current function will organize only the selected WA into `sub_fold` 
% structures, for the purpose of training KRR on only random WA subjects.

ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_ls);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, '/behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

mkdir(fullfile(split_dir, 'randWA'))
for b = 1:nbhvr
    fprintf('#%d behavior: %s\n', b, bhvr_nm{b})
    fn = fullfile(split_dir, 'allAA_randWA', [bhvr_nm{b} split_fstem '.mat']);
    AAWA = load(fn);

    randWA = [];
    for f = 1:length(AAWA.sub_fold)
        randWA = [randWA; AAWA.sub_fold(f).selWA];
    end
    randWA = unique(randWA);
    randWA = intersect(subjects, randWA, 'stable');
    length(randWA)

    CBIG_cell2text(randWA, fullfile(split_dir, 'randWA', ['subj_' bhvr_nm{b} split_fstem '.txt']))
    for f = 1:length(AAWA.sub_fold)
        sub_fold(f).subject_list = AAWA.sub_fold(f).selWA;
        [~, idx] = intersect(randWA, AAWA.sub_fold(f).selWA, 'stable');
        sub_fold(f).fold_index = zeros(length(randWA), 1);
        sub_fold(f).fold_index(idx) = 1;
    end
    save(fullfile(split_dir, 'randWA', [bhvr_nm{b} split_fstem '.mat']), 'sub_fold')
    clear sub_fold
end

end

