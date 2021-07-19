function ABCD_to_HCP_KRR_acc_AAvsWA(KRRoutdir, metric, outmat, ABCD_bhvr_ls, HCP_bhvr_ls)

% ABCD_to_HCP_KRR_acc_AAvsWA(KRRoutdir, metric, outmat, ABCD_bhvr_ls, HCP_bhvr_ls)
%
% 

[ABCD_bhvrs, nbhvr] = CBIG_text2cell(ABCD_bhvr_ls);
HCP_bhvrs = CBIG_text2cell(HCP_bhvr_ls);
assert(length(HCP_bhvrs) == nbhvr, 'Number of behaviors in ABCD_bhvr_ls and HCP_bhvr_ls differ!')

for b = 1:nbhvr
    pred_mat = fullfile(KRRoutdir, ABCD_bhvrs{b}, 'prediction.mat');
    pred = load(pred_mat);
    Nfolds = length(pred.acc);
    N_HCP_iters = size(pred.y_p_AA, 2);

    switch metric
    case 'predictive_COD'
        for f = 1:Nfolds
            for seed = 1:N_HCP_iters
                [pCOD_AA(b,f, seed), pCOD_WA(b,f, seed), pCOD_AAWA(b,f, seed), ss_res_AA(b,f, seed), ...
                    ss_res_WA(b,f, seed), ss_res_AAWA(b,f, seed), ss_total(b,f, seed)] = ABCD_pCOD_2groups...
                    (pred.y_p_AA{f,seed}, pred.y_p_WA{f, seed}, pred.y_t_AA{f, seed}, pred.y_t_WA{f, seed}, pred.y_train{f}, []);
            end
        end
    case 'corr'
        corr_AA(b,:,:) = pred.acc_AA;
        corr_WA(b,:,:) = pred.acc_WA;
    otherwise
        error('Unknown metric: %s', metric)
    end
end

outdir = fileparts(outmat);
if(~exist(outdir, 'dir')); mkdir(outdir); end
switch metric
case 'predictive_COD'
    pCOD_AA = mean(pCOD_AA, 3);  pCOD_WA = mean(pCOD_WA, 3);  pCOD_AAWA = mean(pCOD_AAWA, 3);
    save(outmat, 'pCOD_AA', 'pCOD_WA', 'pCOD_AAWA', 'ss_res_AA', 'ss_res_WA', 'ss_res_AAWA', 'ss_total')
case 'corr'
    corr_AA = mean(corr_AA, 3);  corr_WA = mean(corr_WA, 3);
    save(outmat, 'corr_AA', 'corr_WA')
otherwise
    error('Unknown metric: %s\n', metric)
end
    
end