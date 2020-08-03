function ABCD_PermTest_AAvsWA( group_diff, bhvr_ls, metric, outmat )

alpha = 0.5;

switch metric
    case 'predictive_COD'
        AA_acc = 'pCOD_AA';
        WA_acc = 'pCOD_WA';
    otherwise
        error('Unknown metric.')
end

grpdif = load(group_diff);
acc_diff = grpdif.(WA_acc) - grpdif.(AA_acc);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt';
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
            idx = round(rand(length(grpdif.AA_test{b, f}), 1));
            
            null_yt_AA = grpdif.AA_test{b,f};
            null_yt_AA(idx==1) = grpdif.WA_test{b,f}(idx==1);
            null_yp_AA = grpdif.AA_pred{b,f};
            null_yp_AA(idx==1) = grpdif.WA_pred{b,f}(idx==1);
            
            null_yt_WA = grpdif.WA_test{b,f};
            null_yt_WA(idx==1) = grpdif.AA_test{b,f}(idx==1);
            null_yp_WA = grpdif.WA_pred{b,f};
            null_yp_WA(idx==1) = grpdif.AA_pred{b,f}(idx==1);
            
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
    p_perm(b) = length(find( null_acc_diff(b,:) > abs(avg_acc_diff(b)) | ...
        null_acc_diff < -abs(avg_acc_diff(b)) )) / nperm;
end
p_perm(nbhvr+1) = length(find( mean(null_acc_diff,1) > abs(mean(avg_acc_diff)) | ...
    mean(null_acc_diff,1) < -abs(mean(avg_acc_diff)) )) / nperm;

H_perm_all = FDR(p_perm(:), alpha);
H_perm = setdiff(H_perm_all, nbhvr+1);
sig_diff_idx = sort(H_perm);
sig_diff_bhvr = bhvr_nm(sig_diff_idx);

save(outmat, 'p_perm', 'H_perm_all', 'H_perm', 'sig_diff_idx', 'sig_diff_bhvr');


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
    % remember to revert the sign if using MSE
    otherwise
        error('Unknown metric')
end

end
