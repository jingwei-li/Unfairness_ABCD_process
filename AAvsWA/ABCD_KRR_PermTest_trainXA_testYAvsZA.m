function ABCD_KRR_PermTest_trainXA_testYAvsZA(XAmodel_dir, YA, YA_yvar, ZA, ZA_yvar, bhvr_ls, metric, outmat)

% ABCD_KRR_PermTest_trainXA_testYAvsZA(XAmodel_dir, YA, YA_yvar, ZA, ZA_yvar, bhvr_ls, metric, outmat)
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


nperm = 1000;
alpha = 0.05;
rng('default'); rng(1, 'twister')

%% load accuracy trained on XA
YAacc = []; ZAacc = [];
for b = 1:nbhvr
    fprintf('#%d behavior: %s \n', b, bhvr_nm{b})
    optYA = load(fullfile(XAmodel_dir, bhvr_nm{b}, ['final_result_' YA '_' bhvr_nm{b} '.mat']));
    optZA = load(fullfile(XAmodel_dir, bhvr_nm{b}, ['final_result_' ZA '_' bhvr_nm{b} '.mat']));
    YAacc = [YAacc optYA.optimal_stats.(metric)];
    ZAacc = [ZAacc optZA.optimal_stats.(metric)];

    for f = 1:length(optYA.optimal_stats.(metric))
        for i = 1:nperm
            idx = round(rand(length(optYA.(['yt_' YA_yvar]){f}), 1));

            null_yt_YA = optYA.(['yt_' YA_yvar]){f};
            null_yt_YA(idx==1) = optZA.(['yt_' ZA_yvar]){f}(idx==1);
            null_yp_YA = optYA.(['yp_' YA_yvar]){f};
            null_yp_YA(idx==1) = optZA.(['yt_' ZA_yvar]){f}(idx==1);
            
            null_yt_ZA = optZA.(['yt_' ZA_yvar]){f};
            null_yt_ZA(idx==1) = optYA.(['yt_' YA_yvar]){f}(idx==1);
            null_yp_ZA = optZA.(['yp_' ZA_yvar]){f};
            null_yp_ZA(idx==1) = optYA.(['yp_' YA_yvar]){f}(idx==1);
            
            [null_acc_YA(b,f,i), null_acc_ZA(b,f,i)] = compute_null(null_yt_YA, null_yt_ZA, ...
                null_yp_YA, null_yp_ZA, metric, optYA.y_train_resid{f});   % optYA.y_train_resid is equal to optZA.y_train_resid
        end
    end
end

acc_diff = ZAacc - YAacc;
null_acc_diff = null_acc_ZA - null_acc_YA;
null_acc_diff = squeeze(mean(null_acc_diff, 2));
avg_acc_diff = mean(acc_diff, 1);

p_perm = nan(nbhvr,1);
for b = 1:nbhvr
    length(find( null_acc_diff(b,:) - mean(null_acc_diff(b,:)) > abs(avg_acc_diff(b)) | ...
        null_acc_diff(b,:) - mean(null_acc_diff(b,:)) < -abs(avg_acc_diff(b)) ))
    p_perm(b) = length(find( null_acc_diff(b,:) - mean(null_acc_diff(b,:)) > abs(avg_acc_diff(b)) | ...
        null_acc_diff(b,:) - mean(null_acc_diff(b,:)) < -abs(avg_acc_diff(b)) )) / nperm;
end
length(find( mean(null_acc_diff,1) - mean(mean(null_acc_diff,1),2) > abs(mean(avg_acc_diff)) | ...
    mean(null_acc_diff,1) - mean(mean(null_acc_diff,1),2) < -abs(mean(avg_acc_diff)) ))
p_perm(nbhvr+1) = length(find( mean(null_acc_diff,1) - mean(mean(null_acc_diff,1),2) > abs(mean(avg_acc_diff)) | ...
    mean(null_acc_diff,1) - mean(mean(null_acc_diff,1),2) < -abs(mean(avg_acc_diff)) )) / nperm

H_perm_all = FDR(p_perm(:), alpha)
H_perm = setdiff(H_perm_all, nbhvr+1);
sig_diff_idx = sort(H_perm);
sig_diff_bhvr = bhvr_nm(sig_diff_idx)

save(outmat, 'p_perm', 'H_perm_all', 'H_perm', 'sig_diff_idx', 'sig_diff_bhvr');
    
end



function [null_acc_AA, null_acc_WA] = compute_null(null_yt_AA, null_yt_WA, null_yp_AA, ...
    null_yp_WA, metric, y_train)

switch metric
    case 'predictive_COD'
        ss_res_AA = sum((null_yt_AA - null_yp_AA).^2, 1) ./ length(null_yt_AA);
        ss_res_WA = sum((null_yt_WA - null_yp_WA).^2, 1) ./ length(null_yt_WA);
        ss_total = sum(([null_yt_AA; null_yt_WA] - mean(y_train)).^2, 1) ./ ...
            length([null_yt_AA; null_yt_WA]);
        
        null_acc_AA = bsxfun(@minus, 1, ss_res_AA ./ ss_total);
        null_acc_WA = bsxfun(@minus, 1, ss_res_WA ./ ss_total);
    case 'corr'
        null_acc_AA = CBIG_corr(null_yp_AA, null_yt_AA);
        null_acc_WA = CBIG_corr(null_yp_WA, null_yt_WA);
    case 'MAE'
        null_acc_AA = mean(abs(null_yp_AA - null_yt_AA));
        null_acc_WA = mean(abs(null_yp_WA - null_yt_WA));
    case 'MAE_norm'
        null_acc_AA = mean(abs(null_yp_AA - null_yt_AA))/std(y_train);
        null_acc_WA = mean(abs(null_yp_WA - null_yt_WA))/std(y_train);
    case 'MSE'
        null_acc_AA = mean((null_yp_AA - null_yt_AA).^2);
        null_acc_WA = mean((null_yp_WA - null_yt_WA).^2);
    case 'MSE_norm'
        null_acc_AA = mean((null_yp_AA - null_yt_AA).^2)/var(y_train);
        null_acc_WA = mean((null_yp_WA - null_yt_WA).^2)/var(y_train);
    otherwise
        error('Unknown metric')
end

end

