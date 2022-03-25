function ABCD_scatter_BBAdiff_vs_AccDiff(acc_diff, perm_fname, metric, learned_BBA, real_BBA, outname)

% ABCD_scatter_BBAdiff_vs_AccDiff(acc_diff, perm_fname, metric, learned_BBA, real_BBA, outname)
%
% 

%% determine variable names based on metric type
switch metric
case 'predictive_COD'
    y_label = 'Mean predictive COD: WA - AA';
    AA_acc = 'pCOD_AA';
    WA_acc = 'pCOD_WA';
case 'corr'
    y_label = 'Mean correlation accuracy: WA - AA';
    AA_acc = 'corr_AA';
    WA_acc = 'corr_WA';
otherwise
    error('Unknown accuracy metric: %s', metric)
end
x_label = 'Similarity (learned vs real BBA): WA - AA';

%% load accuracy difference
grpdif = load(acc_diff);
load(perm_fname)

%% deal with model-learned & real brain-behavior associations
load(learned_BBA, 'avg_learned_cov')
realBBA = load(real_BBA, 'avg_cov_testAA', 'avg_cov_testWA');
% check if #behaviors consistent
assert(size(grpdif.(AA_acc), 1) == size(avg_learned_cov,3), ...
    '#behaviors in learned_BBA differred from that in acc_diff')
assert(size(grpdif.(AA_acc), 1) == size(realBBA.avg_cov_testAA,3), ...
    '#behaviors in real_BBA differred from that in acc_diff')
% extract lower-triangular endices
tril_ind = find(tril(ones(size(avg_learned_cov,1)), -1)==1);
learnedBBA = reshape(avg_learned_cov, size(avg_learned_cov,1)*size(avg_learned_cov,2), ...
    size(avg_learned_cov, 3));
realBBA_AA = reshape(realBBA.avg_cov_testAA, size(realBBA.avg_cov_testAA,1) * ...
    size(realBBA.avg_cov_testAA,2), size(realBBA.avg_cov_testAA,3));
realBBA_WA = reshape(realBBA.avg_cov_testWA, size(realBBA.avg_cov_testWA,1) * ...
    size(realBBA.avg_cov_testWA,2), size(realBBA.avg_cov_testWA,3));
learnedBBA = learnedBBA(tril_ind, :);
realBBA_AA = realBBA_AA(tril_ind, :);
realBBA_WA = realBBA_WA(tril_ind, :);

% compute difference in similarity(learned vs real BBA)
simi_AA = zeros(size(learnedBBA,2), 1);  simi_WA = simi_AA;
for b = 1:size(learnedBBA,2)
    simi_AA(b) = CBIG_corr(learnedBBA(:,b), realBBA_AA(:,b));
    simi_WA(b) = CBIG_corr(learnedBBA(:,b), realBBA_WA(:,b));
end
%Xdata = simi_WA(sig_diff_idx) - simi_AA(sig_diff_idx);
%Ydata = mean(grpdif.(WA_acc)(sig_diff_idx,:) - ...
%    grpdif.(AA_acc)(sig_diff_idx,:), 2);
Xdata = simi_WA - simi_AA;
Ydata = mean(grpdif.(WA_acc) - grpdif.(AA_acc), 2);

%% scatter plot
sz = 25;
colormat = [211 94 96] ./ 255;
f = figure;
scatter(Xdata, Ydata, sz, colormat, 'filled')
hold on
xli = get(gca, 'xlim');
xpoints = xli(1):((xli(2)-xli(1))/5):xli(2);

p = polyfit(Xdata, Ydata, 1);
r = polyval(p, xpoints);
plot(xpoints, r, 'k', 'LineWidth', 2)
hold off

[rho, pval] = corr(Xdata, Ydata);
title(sprintf('X- vs Y-axes correlation: %.3f, p value: %e', rho, pval))

xlabel(x_label, 'fontsize', 12)
ylabel(y_label, 'fontsize', 12);

export_fig(outname, '-png', '-nofontswap', '-a1');
set(gcf, 'color', 'w')
hgexport(f, outname)
close
    
end