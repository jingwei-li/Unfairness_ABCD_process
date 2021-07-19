function HCP_to_ABCD_PermTest_AAvsWA(perm_out, group_diff, KRRoutdir, max_HCP_seed, metric, HCP_bhvr_ls, ABCD_bhvr_ls)

% HCP_to_ABCD_PermTest_AAvsWA(perm_out, group_diff, KRRoutdir, max_HCP_seed, metric, HCP_bhvr_ls, ABCD_bhvr_ls)
%
% 

alpha = 0.05;
nperm = 1000;
grpdif = load(group_diff);
switch metric
case 'predictive_COD'
    AA_acc = 'pCOD_AA';
    WA_acc = 'pCOD_WA';
case 'corr'
    AA_acc = 'corr_AA';
    WA_acc = 'corr_WA';
otherwise
    error('Unknown metric: %s', metric)
end

acc_diff = grpdif.(WA_acc) - grpdif.(AA_acc);

[HCP_bhvr, nbhvr] = CBIG_text2cell(HCP_bhvr_ls);
ABCD_bhvr = CBIG_text2cell(ABCD_bhvr_ls);
if(length(ABCD_bhvr) ~= nbhvr)
    error('Numbers of behavioral measures in HCP_bhvr_ls and ABCD_bhvr_ls are different.')
end

rng('default'); rng(1, 'twister')
null_acc_AA = nan(40, nbhvr, nperm);
null_acc_WA = nan(40, nbhvr, nperm);
for b = 1:nbhvr
    c = 1;
    for seed = 1:max_HCP_seed
        pred_mat = fullfile(KRRoutdir, ['randseed_' num2str(seed)], HCP_bhvr{b}, 'prediction_AAWA.mat');
        if(~exist(pred_mat, 'file')); continue; end

        pred = load(pred_mat);
        Nfolds = length(pred.y_t_AA);

        for i = 1:nperm
            % indices of subjects to exchange ethnic/racial group labels
            idx = round(rand(length(pred.y_t_AA{1}{1}), 1));
            % randomly shuffle the ordering of WA subjects
            shuffle = datasample(1:length(pred.y_t_WA{1}{1}), length(pred.y_t_WA{1}{1}), 'replace', false);
            
            curr_null_acc_AA = zeros(Nfolds, 1);
            curr_null_acc_WA = zeros(Nfolds, 1);
            for f = 1:Nfolds
                null_yt_AA = pred.y_t_AA{f}{1};
                null_yt_AA(idx==1) = pred.y_t_WA{f}{1}(shuffle(idx==1));
                null_yp_AA = pred.y_p_AA{f}{1};
                null_yp_AA(idx==1) = pred.y_p_WA{f}{1}(shuffle(idx==1));

                null_yt_WA = pred.y_t_WA{f}{1};
                null_yt_WA(shuffle(idx==1)) = pred.y_t_AA{f}{1}(idx==1);
                null_yp_WA = pred.y_p_WA{f}{1};
                null_yp_WA(shuffle(idx==1)) = pred.y_p_AA{f}{1}(idx==1);

                [curr_null_acc_AA(f), curr_null_acc_WA(f)] = compute_null(null_yt_AA, null_yt_WA, ...
                    null_yp_AA, null_yp_WA, metric, pred.y_train{f}, []);
            end
            null_acc_AA(c,b,i) = mean(curr_null_acc_AA);
            null_acc_WA(c,b,i) = mean(curr_null_acc_WA);
        end
    end
end

null_acc_diff = null_acc_WA - null_acc_AA;
null_acc_diff = squeeze(nanmean(null_acc_diff, 1));
avg_acc_diff = nanmean(acc_diff);

p_perm = nan(nbhvr,1);
for b = 1:nbhvr
    p_perm(b) = length(find( null_acc_diff(b,:) - mean(null_acc_diff(b,:)) > abs(avg_acc_diff(b)) | ...
        null_acc_diff(b,:) - mean(null_acc_diff(b,:)) < -abs(avg_acc_diff(b)) )) / nperm;
end
p_perm(nbhvr+1) = length(find( mean(null_acc_diff,1) - mean(mean(null_acc_diff,1),2) > abs(mean(avg_acc_diff)) | ...
    mean(null_acc_diff,1) - mean(mean(null_acc_diff,1),2) < -abs(mean(avg_acc_diff)) )) / nperm;

H_perm_all = FDR(p_perm(:), alpha);
H_perm = setdiff(H_perm_all, nbhvr+1);
sig_diff_idx = sort(H_perm);
sig_diff_HCP_bhvr = HCP_bhvr(sig_diff_idx);
sig_diff_ABCD_bhvr = ABCD_bhvr(sig_diff_idx);

save(perm_out, 'p_perm', 'H_perm_all', 'H_perm', 'sig_diff_idx', 'sig_diff_HCP_bhvr', ...
    'sig_diff_ABCD_bhvr', 'null_acc_diff', 'avg_acc_diff')
    
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
        null_acc_AA = mean(abs(null_yp_AA - null_yt_AA))/std(AA_train);
        null_acc_WA = mean(abs(null_yp_WA - null_yt_WA))/std(WA_train);
    case 'MSE'
        null_acc_AA = mean((null_yp_AA - null_yt_AA).^2);
        null_acc_WA = mean((null_yp_WA - null_yt_WA).^2);
    case 'MSE_norm'
        null_acc_AA = mean((null_yp_AA - null_yt_AA).^2)/var(AA_train);
        null_acc_WA = mean((null_yp_WA - null_yt_WA).^2)/var(WA_train);
    otherwise
        error('Unknown metric')
end

end
