function  ABCD_whisker_trainXA_testYAvsZA(XAmodel_dir, YA, ZA, bhvr_ls, colloq_ls, metric, tit, outdir, outstem, perm_fname)

% ABCD_whisker_trainXA_testYAvsZA(XAmodel_dir, YA, ZA, bhvr_ls, colloq_ls, metric, tit, outdir, outstem)
%
% Inputs:
%   - XAmodel_dir
%     Directory to KRR models trained on XA subjects (full path).
%   - YA
%     String, the label of YA group, e.g. matchedAA. The test results of YA should be stored in
%     fullfile(XAmodel_dir, <behavior>, ['final_result_' YA '_' <behavior> '.mat'])
%   - ZA 
%     String, the label of ZA group, e.g. matchedWA. The test results of ZA should be stored in 
%     fullfile(XAmodel_dir, <behavior>, ['final_result_' ZA '_' <behavior> '.mat'])
%   - bhvr_ls (optional)
%     Behavior list (full path). Default:
%     '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
%   - colloq_ls (optional)
%     List of collquial name of behaviors (full path). Default:
%     '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt'
%   - metric
%     String, accuracy metric. Choose from 'predictive_COD', 'COD', 'corr', 'MSE' (capital sensitive).
%   - tit
%     Title of the plot.
%   - outdir
%     Output directory (full path).
%   - outstem
%     String, stem of output figure name.
%   - perm_fname
%     A .mat file contains the permutation testing results of AA-WA accuracy difference (full path).
%     It is the result of `ABCD_KRR_PermTest_trainXA_testYAvsZA.m`.

switch metric
    case 'predictive_COD'
        y_label = 'Cross-validated predictive COD';
        y_label_avg = 'Mean cross-validated predictive COD';
    case 'COD'
        y_label = 'Cross-validated COD';
        y_label_avg = 'Mean cross-validated COD';
    case 'corr'
        y_label = 'Cross-validated Pearson''s r';
        y_label_avg = 'Mean cross-validated Pearson''s r';
    case 'MSE'
        y_label = 'Cross-validated MSE';
        y_label_avg = 'Mean cross-validated MSE';
    otherwise
        error('Unknown metric.')
end

if(strfind(YA, 'AA') && strfind(ZA, 'WA'))
    colormat = [114 147 203; 132 186 91; 211 94 96]./255;
elseif(strfind(YA, 'WA') && strfind(ZA, 'AA'))
    colormat = [132 186 91; 114 147 203; 211 94 96]./255;
end
legend1 = YA; legend1(1) = upper(legend1(1));
legend2 = ZA; legend2(1) = upper(legend2(1));
legends = {legend1, legend2, 'Difference'};

%% parse input arguments
ls_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
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

if(exist('perm_fname', 'var') && ~isempty(perm_fname))
    load(perm_fname);
end

%% load accuracy trained on XA
YAacc = []; ZAacc = [];
for b = 1:nbhvr
    optYA = load(fullfile(XAmodel_dir, bhvr_nm{b}, ['final_result_' YA '_' bhvr_nm{b} '.mat']));
    optZA = load(fullfile(XAmodel_dir, bhvr_nm{b}, ['final_result_' ZA '_' bhvr_nm{b} '.mat']));
    if(strcmpi(metric, 'predictive_COD'))
        YAacc = [YAacc optYA.optimal_stats.pCOD];
        ZAacc = [ZAacc optZA.optimal_stats.pCOD];
    else
        YAacc = [YAacc optYA.optimal_stats.(metric)];
        ZAacc = [ZAacc optZA.optimal_stats.(metric)];
    end
end
CV = size(YAacc, 1);

acc_diff = ZAacc - YAacc;
alldata = cat(1, reshape(YAacc, 1, CV, nbhvr), reshape(ZAacc, 1, CV, nbhvr), ...
    reshape(acc_diff, 1, CV, nbhvr));
[~, idx] = sort(mean(acc_diff, 1), 'descend');
data_sort = alldata(:,:,idx);
bhvr_nm_sort = bhvr_nm(idx);
colloq_nm_sort = colloq_nm(idx);

%% if using MSE or MAE metrics: divide all values by the mean of YA MSE in the whole-population model
if(strcmpi(metric, 'MSE') || strcmpi(metric, 'MAE'))
	data_sort_plot = bsxfun(@times, data_sort, 1./mean(data_sort(1,:,:), 2));
else
	data_sort_plot = data_sort;
end

% get indices of significant behaviors
if(exist('perm_fname', 'var') && ~isempty(perm_fname))
    [~, IA, IB] = intersect(idx, sig_diff_idx, 'stable');
else
    IA = [];
end

%% plot for each behavior
ABCD_whisker_2grp_indiv(data_sort_plot, colormat, y_label, legends, ...
    tit, colloq_nm_sort, IA, outdir, outstem, metric);
    
%% plot average across behaviors
avg_data = squeeze(mean(data_sort_plot,2))';
ABCD_whisker_2grp_avg(avg_data, colormat, y_label_avg, legends, outdir, outstem)
    
end