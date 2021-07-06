function ABCD_LRR_predictable_behavior_step1(model_dir, Nperm, subj_ls, csvname)

% ABCD_LRR_predictable_behavior_step1(model_dir, Nperm, subj_ls, csvname)
%
% Wrapper to generate multi-level block permutations.
%
% Inputs:
%   - model_dir
%     Directory of LRR results.
%   - Nperm
%     Number of permutations.
%   - subj_ls  (optional)
%     Subject list.
%     Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt
%   - csvname  (optional)
%     Name of the csv file containing site information.
%     Default: /home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/phenotypes_pass_rs.txt
%
% Output:
%   A file named as [model_dir, 'perm', 'Pset.mat'] will be saved, storing the generated permutations.
%
% Author: Jingwei Li

%% setup default parameters
ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end
alpha_FDR = 0.05;

%% read text files
[subjects, nsub] = CBIG_text2cell(subj_ls);
d = readtable(csvname);
subj_hdr = 'subjectkey';
site_hdr = 'site';
fam_hdr = 'family_id';
grp_hdr = 'rel_grp_id';

[~, ~, idx] = intersect(subjects, d.(subj_hdr), 'stable');
site = d.(site_hdr)(idx);
fam = d.(fam_hdr)(idx);
grp = d.(grp_hdr)(idx);

mkdir(fullfile(model_dir, 'perm'))

%% generate exchangeable blocks and permutations
idx_perm_out = fullfile(model_dir, 'perm', 'Pset.mat');
if(~exist(idx_perm_out, 'file'))
    [Pset, B] = ABCD_multilevel_block_perm(site, fam, grp, Nperm+1);
    N_uq_perm = size(unique(Pset', 'rows'), 1);
    if (N_uq_perm-1 < Nperm)
        warning('Number of required permutations is higher than maximal possible permutations.')
    end
    save(idx_perm_out, 'Pset', 'B')
end

