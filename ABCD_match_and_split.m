function ABCD_match_and_split(csvname, subj_list, race_hdr, site_hdr, fam_hdr, cfds_ls, bhvr_ls, niter, ss_cost_ub, outdir, outstem)

% ABCD_match_and_split(csvname, subj_list, race_hdr, site_hdr, fam_hdr, cfds_ls, bhvr_ls, niter, outdir, outstem)
%
% Inputs:
% - csvname
%   Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/phenotypes_pass_rs.txt'
%
% - subj_list
%   Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt'
%
% - race_hdr
%   Default: 'race'
%
% - site_hdr
%   Default: 'site'
%
% - fam_hdr
%   Default: 'family_id'
%
% - cfds_ls
%   Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/confounds_list.txt'
%
% - bhvr_ls
%   Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
%
% - niter   default: 100
%
% - ss_cost_ub   default: 2.45
%
% Example:
% ABCD_match_and_split([], [], [], [], [], [], [], [], [], ...
%    '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/mat/matchANDsplit/20200719', ...
%    '_pass_rs_pass_pheno')

ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';

if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_list);
subj_hdr = 'subjectkey';

if(~exist('race_hdr', 'var') || isempty(race_hdr))
    race_hdr = 'race';
end

if(~exist('site_hdr', 'var') || isempty(site_hdr))
    site_hdr = 'site';
end
if(~exist('fam_hdr', 'var') || isempty(fam_hdr))
    fam_hdr = 'family_id';
end

if(~exist('cfds_ls', 'var') || isempty(cfds_ls))
    cfds_ls = fullfile(ls_dir, 'confounds_list.txt');
end
cfds_nm = CBIG_text2cell(cfds_ls);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
bhvr_nm = CBIG_text2cell(bhvr_ls);

if(~exist('niter', 'var') || isempty(niter))
    niter = 100;
end

if(~exist('ss_cost_ub', 'var') || isempty(ss_cost_ub))
    ss_cost_ub = 2.45;
end

Nfolds = 10;

%% parse table
d = readtable(csvname);
row_idx = zeros(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects{s});
    if(all(tmp_idx==0))
        error('Subject %s does not exist in %s.', subjects{s}, csvname);
    else
        row_idx(s) = find(tmp_idx==1);
    end
end
d = d(row_idx,:);

race = d.(race_hdr);
site = d.(site_hdr);
for s = 1:nsub
    site{s} = str2num(site{s}(5:end));
end
site = cell2mat(site);

fam_id = d.(fam_hdr);

cfds = [];
n = strcmpi(cfds_nm, 'age');
if(any(n))
    age = d.(cfds_nm{n==1});
    cfds = [cfds age];
end

n = strcmpi(cfds_nm, 'sex');
if(any(n))
    sex = strcmp(d.(cfds_nm{n==1}), 'F');
    cfds = [cfds sex];
end

possib_cfds = {'FD', 'DVARS', 'ICV', 'peduc_avg'};
for c = 1:length(possib_cfds)
    n = strcmpi(cfds_nm, possib_cfds{c});
    if(any(n))
        col = d.(cfds_nm{n==1});
        cfds = [cfds col];
    end
end

cfds_zn = zscore(cfds, 0, 1);

bhvr = zeros(nsub, length(bhvr_nm));
for b = 1:length(bhvr_nm)
    bhvr(:,b) = d.(bhvr_nm{b});
end
bhvr_zn = zscore(bhvr, 0, 1);

%% select AA that can find matched AA, return selected AA and WA
[selAA, selWA, sel_mAA, sel_mWA] = ABCD_match_WAtoAA_within_site(subjects, race, site, bhvr_zn, cfds_zn, niter, ss_cost_ub);

%% statistical testing of measure difference between selected AA and WA
[p_tt, H_tt_FDR, FDR_th] = ABCD_stats_MatchDiff_AAvsWA(sel_mAA, sel_mWA);

%% split selected AA into K folds, as balance as possible
fold_AA = cell(length(bhvr_nm), 1);
for b = 1:length(bhvr_nm)
    curr_selAA = cat(2, selAA{b,:});
    [curr_selAA, curr_AAidx] = intersect(subjects, curr_selAA, 'stable');
    curr_AAsite = site(curr_AAidx);
    curr_AA_fam_id = fam_id(curr_AAidx);
    [fold_AA{b}, fold_sites{b}] = ABCD_split_folds(curr_selAA, curr_AAsite, curr_AA_fam_id, Nfolds);
end

%% based on selected AA folds, split other subjects
fold_subj = cell(length(bhvr_nm), 1);
fold_WA = cell(length(bhvr_nm), 1);
for b = 1:length(bhvr_nm)
    [fold_subj{b}, fold_WA{b}] = ABCD_split_unselected(subjects, site, selAA(b,:), selWA(b,:), fold_AA{b});
end


%% output
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
oname = fullfile(outdir, ['sel_AAWA' outstem '.mat']);
save(oname, 'selAA', 'selWA', 'sel_mAA', 'sel_mWA', 'p_tt', 'H_tt_FDR', 'FDR_th', 'bhvr_nm')

C = combnk(1:Nfolds, 3);
for b = 1:length(bhvr_nm)
    for i = 1:size(C,1)
        sub_fold(i).subject_list = cat(1, fold_subj{b}{C(i,:)});
        sub_fold(i).selAA = intersect(subjects, cat(1, fold_AA{b}{C(i,:)}), 'stable');
        sub_fold(i).selWA = intersect(subjects, cat(1, fold_WA{b}{C(i,:)}), 'stable');
        [sub_fold(i).subject_list, fold_idx] = intersect(subjects, sub_fold(i).subject_list, 'stable');
        sub_fold(i).fold_index = zeros(length(subjects), 1);
        sub_fold(i).fold_index(fold_idx) = 1;
    end
    
    oname = fullfile(outdir, ['sub_fold' outstem '_' bhvr_nm{b} '.mat']);
    save(oname, 'sub_fold');
    
    clear sub_fold
end



end

