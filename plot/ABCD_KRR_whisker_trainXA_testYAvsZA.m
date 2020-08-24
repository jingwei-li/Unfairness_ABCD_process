function  ABCD_whisker_trainXA_testYAvsZA(XAmodel_dir, YA, ZA, bhvr_ls, colloq_ls, metric, tit, outdir, outstem)

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
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
%   - colloq_ls (optional)
%     List of collquial name of behaviors (full path). Default:
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt'
%   - metric
%     String, accuracy metric. Choose from 'predictive_COD', 'COD', 'corr', 'MSE' (capital sensitive).
%   - tit
%     Title of the plot.
%   - outdir
%     Output directory (full path).
%   - outstem
%     String, stem of output figure name.

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

%% load accuracy trained on XA
YAacc = []; ZAacc = [];
for b = 1:nbhvr
    optYA = load(fullfile(XAmodel_dir, bhvr_nm{b}, ['final_result_' YA '_' bhvr_nm{b} '.mat']));
    optZA = load(fullfile(XAmodel_dir, bhvr_nm{b}, ['final_result_' ZA '_' bhvr_nm{b} '.mat']));
    YAacc = [YAacc optYA.optimal_stats.(metric)];
    ZAacc = [ZAacc optZA.optimal_stats.(metric)];
end
CV = size(YAacc, 1);

acc_diff = ZAacc - YAacc;
size(YAacc)
alldata = cat(1, reshape(YAacc, 1, CV, nbhvr), reshape(ZAacc, 1, CV, nbhvr), ...
    reshape(acc_diff, 1, CV, nbhvr));
[~, idx] = sort(mean(acc_diff, 1), 'descend');
data_sort = alldata(:,:,idx);
bhvr_nm_sort = bhvr_nm(idx);
colloq_nm_sort = colloq_nm(idx);

[mean(YAacc(:,idx), 1); mean(ZAacc(:,idx), 1); mean(acc_diff(:,idx), 1)]'
[max(YAacc(:,idx), [], 1); max(ZAacc(:,idx), [], 1); max(acc_diff(:,idx), [], 1)]'
[min(YAacc(:,idx), [], 1); min(ZAacc(:,idx), [], 1); min(acc_diff(:,idx), [], 1)]'

%% plot for each behavior
ABCD_whisker_2grp_indiv(data_sort, colormat, y_label, legends, ...
    tit, colloq_nm_sort, [], outdir, outstem, metric);
    
%% plot average across behaviors
avg_data = squeeze(mean(data_sort,2))';
ABCD_whisker_2grp_avg(avg_data, colormat, y_label_avg, legends, outdir, outstem)
    
end