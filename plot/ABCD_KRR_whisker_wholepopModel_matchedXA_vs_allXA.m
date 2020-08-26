function ABCD_KRR_whisker_wholepopModel_matchedXA_vs_allXA(acc_matched, model_dir, XA, metric, outdir, bhvr_ls, colloq_ls)

% ABCD_KRR_whisker_wholepopModel_matchedXA_vs_allXA(acc_matched, model_dir, XAstem, metric, outdir, bhvr_ls, colloq_ls)
%
% Long description

%% default value for optional argument
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
colloq_nm = CBIG_text2cell(colloq_ls);
assert(length(colloq_nm)==nbhvr, 'Numbers of behaviors in colloq_ls and bhvr_ls are different.')

%% load accuracy of matched sub-sample
acc_matched = load(acc_matched);

switch metric
case 'predictive_COD'
    acc_prefix = 'pCOD';
otherwise
    acc_prefix = metric;
end
acc_field = [acc_prefix '_' XA];
acc_matched = acc_matched.(acc_field);
assert(size(acc_matched, 1) == nbhvr, 'Numbers of behaviors in bhvr_ls and acc_matched are different')

%% load accuracy of all subjects in sub-population
acc_all = nan(size(acc_matched));
for b = 1:nbhvr
    opt = load(fullfile(model_dir, ['final_result_all' XA '_' bhvr_nm{b} '.mat']));
    acc_all(b,:) = opt.optimal_stats.(metric);
end

%% setup figure utilities
if(strcmp(XA, 'AA'))
	colormat = [114 147 203; 137 216 248; 211 94 96]./255;
elseif (strcmp(XA, 'WA'))
	colormat = [132 186 91; 177 231 179; 211 94 96]./255;
end

switch metric
case 'predictive_COD'
    yl_appendix = 'predictive COD';
case 'corr'
    yl_appendix = 'Pearson''s r';
case {'MSE', 'MAE'}
    yl_appendix = metric;
case {'MSE_norm', 'MAE_norm'}
    yl_appendix = ['normalized ' metric];
otherwise
    error('Unknown metric: %s', metric)
end
y_label = ['Cross-validated ' yl_appendix];
y_label_avg = ['Mean cross-validated ' yl_appendix];

legends = {['Matched ' XA], ['All ' XA], 'Difference'};
tit = sprintf('Whole-population model, accuracy: matched %s vs all %s', XA, XA);
mkdir(outdir)

%% plot for individual behaviors
data = cat(3, acc_matched, acc_all, acc_all - acc_matched);
data = permute(data, [3 2 1]);
ABCD_whisker_2grp_indiv(data, colormat, y_label, legends, ...
    tit, colloq_nm, [], outdir, metric, metric)
    
%% plot for average
avg_data = mean(data, 3)';
ABCD_whisker_2grp_avg(avg_data, colormat, y_label_avg, legends, outdir, metric)
    
end