function HCP_to_ABCD_KRR_acc_AAvsWA(KRRoutdir, metric, outmat, max_HCP_seed, HCP_bhvr_ls, ABCD_bhvr_ls)

% HCP_to_ABCD_KRR_acc_AAvsWA(KRRoutdir, metric, outmat, HCP_bhvr_ls, ABCD_bhvr_ls)
%
% 

[HCP_bhvr, nbhvr] = CBIG_text2cell(HCP_bhvr_ls);
ABCD_bhvr = CBIG_text2cell(ABCD_bhvr_ls);
if(length(ABCD_bhvr) ~= nbhvr)
    error('Numbers of behavioral measures in HCP_bhvr_ls and ABCD_bhvr_ls are different.')
end

for b = 1:nbhvr
    c = 1;
    for seed = 1:max_HCP_seed
        pred_mat = fullfile(KRRoutdir, ['randseed_' num2str(seed)], HCP_bhvr{b}, 'prediction_AAWA.mat');
        if(~exist(pred_mat, 'file')); continue; end

        pred = load(pred_mat);
        switch metric
        case 'predictive_COD'
            for f = 1:length(pred.y_p_AA)
                [pCOD_AA(c,b,f), pCOD_WA(c,b,f), pCOD_AAWA(c,b,f), ss_res_AA(c,b,f), ...
                    ss_res_WA(c,b,f), ss_res_AAWA(c,b,f), ss_total(c,b,f)] = ABCD_pCOD_2groups...
                    (pred.y_p_AA{f}{1}, pred.y_p_WA{f}{1}, pred.y_t_AA{f}{1}, pred.y_t_WA{f}{1}, pred.y_train{f}, []);
            end
        case 'corr'
            corr_AA(c,b,:) = pred.acc_AA;
            corr_WA(c,b,:) = pred.acc_WA;
        otherwise
            error('Unknown metric: %s', metric)
        end
        
        c = c + 1;
    end
end

switch metric
case 'predictive_COD'
    pCOD_AA = mean(pCOD_AA, 3);  pCOD_WA = mean(pCOD_WA, 3);  pCOD_AAWA = mean(pCOD_AAWA, 3);
    ss_res_AA = mean(ss_res_AA, 3);  ss_res_WA = mean(ss_res_WA, 3); ss_res_AAWA = mean(ss_res_AAWA, 3);
    ss_total = mean(ss_total, 3);

    save(outmat, 'pCOD_AA', 'pCOD_WA', 'pCOD_AAWA', 'ss_res_AA', 'ss_res_WA', 'ss_res_AAWA', 'ss_total')
case 'corr'
    corr_AA = mean(corr_AA, 3);  corr_WA = mean(corr_WA, 3);
    save(outmat, 'corr_AA', 'corr_WA')
end
    
end