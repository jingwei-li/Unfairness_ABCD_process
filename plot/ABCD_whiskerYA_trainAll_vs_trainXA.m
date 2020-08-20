function ABCD_whiskerYA_trainAll_vs_trainXA(YA, bhvr_ls, colloq_ls, acc_trainAll, ...
	model_trainXA, metric, tit, outdir, outstem)

% ABCD_whiskerYA_trainAll_vs_trainXA()
%
% Long description

switch metric
	case 'predictive_COD'
		y_label = 'Cross-validated predictive COD';
		y_label_avg = 'Mean cross-validated predictive COD';
		YA_acc = ['pCOD_' YA];
	case 'corr'
		y_label = 'Cross-validated Pearson''s r';
		y_label_avg = 'Mean cross-validated Pearson''s r';
		YA_acc = ['corr_' YA];
	otherwise
		error('Unknown metric.')
end

if(strcmp(YA, 'AA'))
	colormat = [114 147 203; 137 216 248; 211 94 96]./255;
elseif (strcmp(YA, 'WA'))
	colormat = [132 186 91; 177 231 179; 211 94 96]./255;
end
legends = {'whole-pop trained', 'AA trained', 'Difference'};

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
trainOnAll = load(acc_trainAll);
trainOnAll = trainOnAll.(YA_acc);
CV = size(trainOnAll, 2);

%% load accuracy trained on AA
trainOnXA = nan(size(trainOnAll));
for b = 1:nbhvr
	opt = load(fullfile(model_trainXA, bhvr_nm{b}, ['final_result_matched' YA '_' bhvr_nm{b} '.mat']));
	trainOnXA(b,:) = opt.optimal_stats.(metric);
end

acc_diff = trainOnXA - trainOnAll;
alldata = cat(1, reshape(trainOnAll', 1, CV, nbhvr), reshape(trainOnXA', 1, CV, nbhvr), ...
	reshape(acc_diff', 1, CV, nbhvr));
[~, idx] = sort(mean(acc_diff, 2), 'descend');
data_sort = alldata(:,:,idx);
bhvr_nm_sort = bhvr_nm(idx);
colloq_nm_sort = colloq_nm(idx);

%% plot for each behavior
ABCD_whisker_2grp_indiv(data_sort, colormat, y_label, legends, ...
	tit, colloq_nm_sort, [], outdir, outstem)
	
end