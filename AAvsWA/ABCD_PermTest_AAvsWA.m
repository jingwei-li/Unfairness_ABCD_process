function ABCD_PermTest_AAvsWA( group_diff, bhvr_ls, metric, outmat )

% ABCD_PermTest_AAvsWA( group_diff, bhvr_ls, metric, outmat )
%
% Input
%   - group_diff
%     A .mat file contains the accuracy of each race group, and the
%     original & predicted scores of each group (full path).
%
%   - bhvr_ls (optional)
%     Behavior list (full path, text file). Use the behaviors which passed
%     predictability criteria.
%     Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
%
%   - metric
%     Accuracy metric. For now, only support 'predictive_COD', 'corr'.
% 
%   - outmat
%     Output filename (full path).
% 

alpha = 0.5;
grpdif = load(group_diff);

switch metric
    case 'predictive_COD'
        AA_acc = 'pCOD_AA';
        WA_acc = 'pCOD_WA';
    case 'corr'
        AA_acc = 'corr_AA';
        WA_acc = 'corr_WA';
        grpdif.AA_train = cell(size(grpdif.corr_AA));
        grpdif.WA_train = cell(size(grpdif.corr_WA));
    case 'MSE'
        AA_acc = 'MSE_AA';
        WA_acc = 'MSE_WA';
        grpdif.AA_train = cell(size(grpdif.MSE_AA));
        grpdif.WA_train = cell(size(grpdif.MSE_WA));
    otherwise
        error('Unknown metric.')
end

acc_diff = grpdif.(WA_acc) - grpdif.(AA_acc);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt';
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

nperm = 1000;
alpha = 0.05;

%% permutation
rng('default'); rng(1, 'twister')
null_acc_AA = nan(nbhvr, size(acc_diff,2), nperm);
null_acc_WA = nan(nbhvr, size(acc_diff,2), nperm);
for b = 1:nbhvr
    for f = 1:size(acc_diff, 2)
        for i = 1:nperm
            % indices of subjects to exchange ethnic/racial group labels
            idx = round(rand(length(grpdif.AA_test{b, f}), 1));
            % randomly shuffle the ordering of WA subjects
            shuffle = datasample(1:length(grpdif.WA_test{b,f}), length(grpdif.WA_test{b,f}), 'replace', false);
            
            null_yt_AA = grpdif.AA_test{b,f};
            null_yt_AA(idx==1) = grpdif.WA_test{b,f}(shuffle(idx==1));
            null_yp_AA = grpdif.AA_pred{b,f};
            null_yp_AA(idx==1) = grpdif.WA_pred{b,f}(shuffle(idx==1));
            
            null_yt_WA = grpdif.WA_test{b,f};
            null_yt_WA(shuffle(idx==1)) = grpdif.AA_test{b,f}(idx==1);
            null_yp_WA = grpdif.WA_pred{b,f};
            null_yp_WA(shuffle(idx==1)) = grpdif.AA_pred{b,f}(idx==1);
            
            [null_acc_AA(b,f,i), null_acc_WA(b,f,i)] = compute_null(null_yt_AA, null_yt_WA, ...
                null_yp_AA, null_yp_WA, metric, grpdif.AA_train{b,f}, grpdif.WA_train{b,f});
        end
    end
end

null_acc_diff = null_acc_WA - null_acc_AA;
null_acc_diff = squeeze(mean(null_acc_diff, 2));
avg_acc_diff = mean(acc_diff, 2);

p_perm = nan(nbhvr,1);
for b = 1:nbhvr
    p_perm(b) = length(find( null_acc_diff(b,:) - mean(null_acc_diff(b,:)) > abs(avg_acc_diff(b)) | ...
        null_acc_diff(b,:) - mean(null_acc_diff(b,:)) < -abs(avg_acc_diff(b)) )) / nperm;
    %p_perm(b) = length(find( null_acc_diff(b,:)  > abs(avg_acc_diff(b)) | ...
    %    null_acc_diff(b,:)  < -abs(avg_acc_diff(b)) )) / nperm;
end
p_perm(nbhvr+1) = length(find( mean(null_acc_diff,1) - mean(mean(null_acc_diff,1),2) > abs(mean(avg_acc_diff)) | ...
    mean(null_acc_diff,1) - mean(mean(null_acc_diff,1),2) < -abs(mean(avg_acc_diff)) )) / nperm;

H_perm_all = FDR(p_perm(:), alpha);
H_perm = setdiff(H_perm_all, nbhvr+1);
sig_diff_idx = sort(H_perm);
sig_diff_bhvr = bhvr_nm(sig_diff_idx);

save(outmat, 'p_perm', 'H_perm_all', 'H_perm', 'sig_diff_idx', 'sig_diff_bhvr', 'null_acc_diff');


end

function [null_acc_AA, null_acc_WA] = compute_null(null_yt_AA, null_yt_WA, null_yp_AA, ...
    null_yp_WA, metric, AA_train, WA_train)

switch metric
    case 'predictive_COD'
        ss_res_AA = sum((null_yt_AA - null_yp_AA).^2, 1) ./ length(null_yt_AA);
        ss_res_WA = sum((null_yt_WA - null_yp_WA).^2, 1) ./ length(null_yt_WA);
        ss_total = sum(([null_yt_AA; null_yt_WA] - mean([AA_train; WA_train])).^2, 1) ./ ...
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

