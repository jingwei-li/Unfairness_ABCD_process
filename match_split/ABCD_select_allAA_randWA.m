function ABCD_select_allAA_randWA(csvname, race_hdr, site_hdr, bhvr_ls, subj_ls, split_dir, split_fstem)

% ABCD_select_allAA_randWA(csvname, race_hdr, site_hdr, bhvr_ls, subj_ls, split_dir, split_fstem)
%
% Given the folds that have already been split, select all AA subjects, 
% and randomly select the same number of WA subjects.

ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end
if(~exist('race_hdr', 'var') || isempty(race_hdr))
    race_hdr = 'race';
end
if(~exist('site_hdr', 'var') || isempty(site_hdr))
    site_hdr = 'site';
end
subj_hdr = 'subjectkey';

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_ls);
if(size(subjects,1) == 1)
    subjects = subjects';
end

%% read table
d = readtable(csvname);
race = d.(race_hdr);  site = d.(site_hdr);
site = cellfun(@(x)x(5:6), site, 'uniformoutput', false);
site = cellfun(@str2num, site);
all_subj = d.(subj_hdr);
[~, ~, idx] = intersect(subjects, all_subj, 'stable');
race = race(idx);  site = site(idx);
uq_site = unique(site);
AA = subjects(race==2);  siteAA = site(race==2);
WA = subjects(race==1);  siteWA = site(race==1);

%% randomly select WA for each site
rng('default'); rng(1000)
AA_st = cell(length(uq_site), 1);  WA_st = AA_st;
for st = 1:length(uq_site)
    AA_st{uq_site(st)} = AA(siteAA == uq_site(st));
    if(isempty(AA_st{uq_site(st)}))
        WA_st{uq_site(st)} = [];
    else
        WA_st{uq_site(st)} = WA(siteWA == uq_site(st));
        if( length(AA_st{uq_site(st)}) < length(WA_st{uq_site(st)}) )
            idx = datasample(1:length(WA_st{uq_site(st)}), ...
                length(AA_st{uq_site(st)}), 'replace', false);
            WA_st{uq_site(st)} = WA_st{uq_site(st)}(sort(idx, 'ascend'));
        elseif( length(AA_st{uq_site(st)}) > length(WA_st{uq_site(st)}) )
            idx = datasample(1:length(AA_st{uq_site(st)}), ...
                length(WA_st{uq_site(st)}), 'replace', false);
            AA_st{uq_site(st)} = AA_st{uq_site(st)}(sort(idx, 'ascend'));
        end
    end
end

%% for each fold, aggregate selected AA and WA base on corresponding sites
mkdir(fullfile(split_dir, 'allAA_randWA'))
for b = 1:nbhvr
    fprintf('#%d behavior: %s\n', b, bhvr_nm{b})
    fname = fullfile(split_dir, ['sub_fold' split_fstem '_' bhvr_nm{b} '.mat']);
    split = load(fname);

    for f = 1:length(split.sub_fold)
        [~, ~, idx] = intersect(split.sub_fold(f).subject_list, all_subj, 'stable');

        site_fold = site(split.sub_fold(f).fold_index == 1);
        uq_site_fold = unique(site_fold);

        currAA = []; currWA = [];
        for st = 1:length(uq_site_fold)
            addAA = intersect(AA_st{uq_site_fold(st)}, split.sub_fold(f).subject_list);
            addWA = intersect(WA_st{uq_site_fold(st)}, split.sub_fold(f).subject_list);
            currAA = [currAA; addAA];
            currWA = [currWA; addWA];
        end

        sub_fold(f).subject_list = split.sub_fold(f).subject_list;
        sub_fold(f).selAA = currAA;
        sub_fold(f).selWA = currWA;
        sub_fold(f).fold_index = split.sub_fold(f).fold_index;
    end
    outname = fullfile(split_dir, 'allAA_randWA', [bhvr_nm{b} split_fstem '.mat']);
    save(outname, 'sub_fold')
    clear sub_fold
end

end