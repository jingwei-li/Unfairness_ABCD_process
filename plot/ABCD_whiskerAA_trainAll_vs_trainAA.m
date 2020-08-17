function ABCD_whiskerAA_trainAll_vs_trainAA(bhvr_ls, colloq_ls, acc_trainAll, model_trainAA, metric)

% ABCD_whiskerAA_trainAll_vs_trainAA()
%
% Long description

switch metric
case 'predictive_COD'
	y_label = 'Cross-validated predictive COD';
	y_label_avg = 'Mean cross-validated predictive COD';
	AA_acc = 'pCOD_AA';
case 'corr'
	y_label = 'Cross-validated Pearson''s r';
	y_label_avg = 'Mean cross-validated Pearson''s r';
	AA_acc = 'corr_AA';
otherwise
	error('Unknown metric.')
end

colormat = [114 147 203； 40 50 149； 211 94 96]./255;
legends = {'whole-pop trained', 'AA trained', 'Difference'};
tit = 'Compare '

%% parse input arguments
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

%% load accuracy trained on whole population
trainAll = load(acc_trainAll);
trainAll = trainAll.(AA_acc);
CV = size(trainAll, 2);

%% load accuracy trained on AA
trainAA = nan(size(trainAll));
for b = 1:nbhvr
	opt = fullfile(model_trainAA, bhvr_nm{b}, ['final_result_' bhvr_nm{b} '.mat']);
	trainAA(b,:) = opt.optimal_stats.(metric);
end

acc_diff = trainAA - trainAll;
alldata = cat(1, reshape(trainAll', 1, CV, nbhvr), reshape(trainAA', 1, CV, nbhvr), ...
	reshape(acc_diff', 1, CV, nbhvr));
[~, idx] = sort(mean(acc_diff, 1), 'descend');
data_sort = alldata(:,:,idx);
bhvr_nm_sort = bhvr_nm(idx);
colloq_nm_sort = colloq_nm(idx);

%% plot for each behavior
ABCD_whisker_2grp_indiv(data_sort, colormat, y_label, legends, ...
	tit, colloq_nm, sigdiff_idx, outdir, outstem)
	
end