function ABCD_whisker_AAvsWA(bhvr_ls, colloq_ls, group_diff, perm_fname, metric, outdir, outstem)

% ABCD_whisker_AAvsWA(bhvr_ls, colloq_ls, group_diff, perm_fname, metric, outdir, outstem)
% 
% Input:
%   - bhvr_ls (optional)
%     Behavior list (full path, text file). Use the behaviors which passed
%     predictability criteria. Default: 
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
% 
%   - colloq_ls (optional)
%     List of behaviors' colloquial names (full path). Use the same
%     behaviors as in bhvr_ls. Default: 
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt'
%
%   - group_diff
%     A .mat file contains the accuracy of each race group, and the
%     original & predicted scores of each group (full path).
%
%   - perm_fname
%     A .mat file contains the permutation testing results of AA-WA
%     accuracy difference (full path).
%
%   - metric
%     Type of accuracy metric. For now only support 'predictive_COD', 
%     'corr'.
%
%   - outdir
%     Output directory (full path).
%
%   - outstem
%     A string to be atttached onto the output filenames, e.g.
%     '_pass_rs_pass_pheno'.
%


%% figure utilities setup
switch metric
    case 'predictive_COD'
        y_label = 'Cross-validated predictive COD';
        y_label_avg = 'Mean cross-validated predictive COD';
        AA_acc = 'pCOD_AA';
        WA_acc = 'pCOD_WA';
    case 'corr'
        y_label = 'Cross-validated Pearson''s r';
        y_label_avg = 'Mean cross-validated Pearson''s r';
        AA_acc = 'corr_AA';
        WA_acc = 'corr_WA';
    otherwise
        error('Unknown metric.')
end

colormat = [114 147 203; 132 186 91; 211 94 96]./255;
legends = {'AA', 'Matched WA', 'Difference'};

%% parse input arguments, collect data
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
colloq_nm = CBIG_text2cell(colloq_ls);
% replace '_' with ' '
for b = 1:nbhvr
    curr_colloq = strsplit(colloq_nm{b}, '_');
    colloq_nm{b} = strjoin(curr_colloq, ' ');
end

grpdif = load(group_diff);
if(exist('perm_fname', 'var') && ~isempty(perm_fname))
    load(perm_fname);
end

CV = size(grpdif.(AA_acc), 2);

acc_diff = grpdif.(WA_acc)' - grpdif.(AA_acc)';
alldata = cat(1, reshape(grpdif.(AA_acc)', 1, CV, nbhvr), reshape(grpdif.(WA_acc)', 1, CV, nbhvr), ...
    reshape(acc_diff, 1, CV, nbhvr));
[~, idx] = sort(mean(acc_diff, 1), 'descend');
data_sort = alldata(:,:,idx);
bhvr_nm_sort = bhvr_nm(idx);
colloq_nm_sort = colloq_nm(idx);

% get indices of significant behaviors
if(exist('perm_fname', 'var') && ~isempty(perm_fname))
    [~, IA, IB] = intersect(idx, sig_diff_idx, 'stable');
else
    IA = [];
end

if(~exist(outdir, 'dir'))
    mkdir(outdir)
end

%% plot for each behavior
ABCD_whisker_2grp_indiv(data_sort, colormat, y_label, legends, ...
	colloq_nm_sort, IA, outdir, outstem)

%% plot the average
avg_data = mean(data_sort, 3)';
ABCD_whisker_2grp_avg(avg_data, colormat, y_label, legends, outdir, outstem)

end

