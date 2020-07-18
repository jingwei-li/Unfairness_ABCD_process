function ABCD_split_and_match(csvname, subj_list, race_hdr, site_hdr, fam_hdr, cfds_ls, bhvr_ls, niter, outdir, outname)

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

if(~exist('cfds_ls', 'var') || isempty(cfds_hdr))
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

%%
[selAA, selWA, sel_mAA, sel_mWA] = ABCD_match_WAtoAA_within_site(subjects, race, site, bhvr_zn, cfds_zn, niter);

for b = 1:length(bhvr_nm)
    curr_selAA = cat(2, selAA{b,:});
    [curr_selAA, curr_AAidx] = intersect(subjects, curr_selAA, 'stable');
    curr_AAsite = site(curr_AAidx);
    curr_AA_fam_id = fam_id(curr_AAidx);
    [fold_AA{b}, fold_sites{b}] = ABCD_split_folds(curr_selAA, curr_AAsite, curr_AA_fam_id, 9);
end





end

