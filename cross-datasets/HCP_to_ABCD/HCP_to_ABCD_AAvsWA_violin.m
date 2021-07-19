function HCP_to_ABCD_AAvsWA_violin(group_diff, perm_stats, metric, fig_out, ABCD_bhvr_ls, ABCD_colloq_ls)

% HCP_to_ABCD_AAvsWA_violin(acc_mat, perm_stats, metric, fig_out, ABCD_bhvr_ls, ABCD_colloq_ls)
%
% 

%% figure utilities setup
switch metric
case 'predictive_COD'
    y_label = 'Cross-validated predictive COD';
    y_label_avg = 'Mean cross-validated predictive COD';
    AA_acc = 'pCOD_AA';
    WA_acc = 'pCOD_WA';
case 'corr'
    y_label = 'Cross-validated Pearson''s r';
    y_label_avg = 'Mean cross-validated Pearson''s r';
    AA_acc = 'corr_AA';
    WA_acc = 'corr_WA';
case 'MSE'
    y_label = 'Cross-validated MSE';
    y_label_avg = 'Mean cross-validated MSE';
    AA_acc = 'MSE_AA';
    WA_acc = 'MSE_WA';
otherwise
    error('Unknown metric.')
end

grpdif = load(group_diff);
stats = load(perm_stats);

colormat = [114 147 203; 132 186 91; 211 94 96]./255;
legends = {'AA', 'Matched WA', 'Difference'};

[bhvr_nm, nbhvr] = CBIG_text2cell(ABCD_bhvr_ls);
colloq_nm = CBIG_text2cell(ABCD_colloq_ls);
assert(length(colloq_nm) == nbhvr, ...
    'Difference lengths of behavioral list and colloquial name list.')
% replace '_' with ' '
for b = 1:nbhvr
    curr_colloq = strsplit(colloq_nm{b}, '_');
    colloq_nm{b} = strjoin(curr_colloq, ' ');
end

nseeds = size(grpdif.(AA_acc), 1);
acc_diff = grpdif.(WA_acc) - grpdif.(AA_acc);
alldata = cat(1, reshape(grpdif.(AA_acc), 1, nseeds, nbhvr), reshape(grpdif.(WA_acc), 1, nseeds, nbhvr), ...
    reshape(acc_diff, 1, nseeds, nbhvr));

%% if using MSE or MAE metrics: divide all values by the mean of AA MSE
if(strcmpi(metric, 'MSE') || strcmpi(metric, 'MAE'))
	data_plot = bsxfun(@times, alldata, 1./mean(alldata(1,:,:), 2));
else
	data_plot = alldata;
end

[outdir, outstem] = fileparts(fig_out);
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
%% plot for each behavior
ABCD_violin_2grp_indiv(data_plot, colormat, y_label, legends, [], ...
	colloq_nm, stats.sig_diff_idx, outdir, outstem, false)
    
end