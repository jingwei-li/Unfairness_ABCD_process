function ABCD_KRR_whisker_models_wholepop_vs_subgrp(wholepop_dir, subgrp_dir, legends, outdir, bhvr_ls, colloq_ls)

% ABCD_KRR_whisker_models_wholepop_vs_subgrp(wholepop_dir, subgrp_dir, legends, outdir, bhvr_ls, coolq_ls)
%
% Long description

ls_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
[colloq_nm] = CBIG_text2cell(colloq_ls);

%% obtain all possible metrics
tmp_opt = load(fullfile(wholepop_dir, ['final_result_' bhvr_nm{1} '.mat']));
metrics = fieldnames(tmp_opt.optimal_stats);
CV = length(tmp_opt.optimal_stats.(metrics{1}));

%% collect accuracies of every metric for each model
for b = 1:nbhvr
    wholepop_opt = load(fullfile(wholepop_dir, ['final_result_' bhvr_nm{b} '.mat']));
    subgrp_opt = load(fullfile(subgrp_dir, bhvr_nm{b}, ['final_result_' bhvr_nm{b} '.mat']));

    for m = 1:length(metrics)
        wholepop_acc.(metrics{m})(:,b) = wholepop_opt.optimal_stats.(metrics{m});
        subgrp_acc.(metrics{m})(:,b) = subgrp_opt.optimal_stats.(metrics{m});
        acc_diff.(metrics{m})(:,b) = subgrp_acc.(metrics{m})(:,b) - wholepop_acc.(metrics{m})(:,b);
    end
end

%% plot for each metric
colormat = [255 255 255; 192 192 192; 128 128 128]./255;
mkdir(outdir)
for m = 1:length(metrics)
    y_label = 'Cross-validated ';
    y_label_avg = 'Mean cross-validated ';
    switch metrics{m}
    case 'corr'
        y_label = [y_label 'Pearson''s r'];
        y_label_avg = [y_label_avg 'Pearson''s r'];
    case 'predictive_COD'
        y_label = [y_label 'predictive COD'];
        y_label_avg = [y_label_avg 'predictive COD'];
    case 'MAE_norm'
        y_label = [y_label 'normalized MAE'];
        y_label_avg = [y_label_avg 'normalized MAE'];
    case 'MSE_norm'
        y_label = [y_label 'normalized MSE'];
        y_label_avg = [y_label_avg 'normalized MSE'];
    otherwise
        y_label = [y_label metrics{m}];
        y_label_avg = [y_label_avg metrics{m}];
    end
    
    data = cat(1, reshape(wholepop_acc.(metrics{m}), 1, CV, nbhvr), ...
        reshape(subgrp_acc.(metrics{m}), 1, CV, nbhvr), reshape(acc_diff.(metrics{m}), 1, CV, nbhvr));
    ABCD_whisker_2grp_indiv(data, colormat, y_label, legends, ...
        [], colloq_nm, [], outdir, metrics{m}, metrics{m})
    
    avg_data = mean(data, 3)';
    ABCD_whisker_2grp_avg(avg_data, colormat, y_label_avg, legends, outdir, metrics{m})
end
    
end