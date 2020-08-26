function ABCD_create_subfold_allAA_randWA(subj_ls, bhvr_ls, split_dir, split_fstem)

% ABCD_create_subfold_allAA_randWA(subj_ls, bhvr_ls, split_dir, split_fstem)
%
% This function is based on the results from `ABCD_select_allAA_randWA.m`.
% `ABCD_select_allAA_randWA.m` has selected all AA and random WA with the 
% same sample size (if #WA < #AA for some site, only part of AA subjects 
% are selected to match the size).
% The current function will organize only the selected AA and WA into `sub_fold` 
% structures, for the purpose of training KRR on only selected AA and WA subjects.

ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_ls);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, '/behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

out_subdir = 'train_allAA_randWA';
mkdir(fullfile(split_dir, out_subdir))
for b = 1:nbhvr
    fprintf('#%d behavior: %s\n', b, bhvr_nm{b})
    fn = fullfile(split_dir, 'allAA_randWA', [bhvr_nm{b} split_fstem '.mat']);
    AAWA = load(fn);

    selected = [];
    for f = 1:length(AAWA.sub_fold)
        selected = [selected; AAWA.sub_fold(f).selAA; AAWA.sub_fold(f).selWA];
    end
    selected = unique(selected);
    selected = intersect(subjects, selected, 'stable');
    length(selected)

    CBIG_cell2text(selected, fullfile(split_dir, out_subdir, ['subj_' bhvr_nm{b} split_fstem '.txt']))
    for f = 1:length(AAWA.sub_fold)
        [sub_fold(f).subject_list, idx] = intersect(selected, [AAWA.sub_fold(f).selAA; AAWA.sub_fold(f).selWA], 'stable');
        sub_fold(f).fold_index = zeros(length(selected), 1);
        sub_fold(f).fold_index(idx) = 1;
    end
    save(fullfile(split_dir, out_subdir, [bhvr_nm{b} split_fstem '.mat']), 'sub_fold')
    clear sub_fold
end
    
end