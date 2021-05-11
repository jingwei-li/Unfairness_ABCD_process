function ABCD_KRR_yt_yp_diff_AAvsWA(colloq_ls, colloq_of_interest, pCOD_diff, outmat)


    %% parse input arguments
    ls_dir = '/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
    if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
        colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
    end
    [colloq_nm, nbhvr] = CBIG_text2cell(colloq_ls);

    [COI] = CBIG_text2cell(colloq_of_interest);
    [~,~,c_idx] = intersect(COI, colloq_nm, 'stable');

    
    grpdif = load(pCOD_diff);
    CV = size(grpdif.AA_test, 2);
    
    if(size(grpdif.AA_test, 1) ~= nbhvr)
        error('Numbers of behaviors in colloq_ls and group_diff are different')
    end
    
    %% scatter plot
    corr_AA = nan(nbhvr, CV); corr_WA = corr_AA;
    shift_AA = nan(nbhvr, CV); shift_WA = shift_AA;
    shift_sq_AA = nan(nbhvr, CV); shift_sq_WA = shift_sq_AA;
    for b = 1:nbhvr
        fprintf('#%d behavior: %s\n', b, colloq_nm{b})
    
        for cv = 1:CV
            corr_AA(b,cv) = CBIG_corr(grpdif.AA_test{b,cv}, grpdif.AA_pred{b,cv});
            corr_WA(b,cv) = CBIG_corr(grpdif.WA_test{b,cv}, grpdif.WA_pred{b,cv});
    
            shift_AA(b,cv) = mean(grpdif.AA_test{b,cv} - grpdif.AA_pred{b,cv});
            shift_WA(b,cv) = mean(grpdif.WA_test{b,cv} - grpdif.WA_pred{b,cv});
    
            shift_sq_AA(b,cv) = shift_AA(b,cv)^2;
            shift_sq_WA(b,cv) = shift_WA(b,cv)^2;
        end
        p(b) = ttest2(shift_sq_AA(b,:), shift_sq_WA(b,:), 'vartype', 'unequal');
    end
    [H, th] = FDR(p(c_idx), 0.05);
    disp(COI(H))
    
    save(outmat, 'p', 'COI', 'H', 'shift_AA', 'shift_WA', 'shift_sq_AA', 'shift_sq_WA')
    end
    
    