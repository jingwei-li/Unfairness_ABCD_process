function HCP_to_ABCD_violin_allsubj(test_dir, metric, outdir, outstem, max_HCP_seed, HCP_bhvr_ls, ABCD_colloq_ls)

% HCP_to_ABCD_general_acc_violin(test_dir, metric, outdir, outstem, max_HCP_seed, HCP_bhvr_ls, ABCD_colloq_ls)
%
% 

addpath(genpath( fullfile(getenv('HOME'), 'storage', 'from_HOME', 'code', 'plotting_functions')))

[HCP_bhvrs, nbhvr] = CBIG_text2cell(HCP_bhvr_ls);
ABCD_colloqs = CBIG_text2cell(ABCD_colloq_ls);
assert(length(ABCD_colloqs) == nbhvr, 'Number of behaviors in HCP_bhvr_ls and ABCD_colloq_ls differ!')

metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
[~, midx] = intersect(metrics, metric, 'stable');
for b = 1:nbhvr
    c = 0;
    for seed = 1:max_HCP_seed
        pred_file = fullfile(test_dir, ['randseed_' num2str(seed)], HCP_bhvrs{b}, 'prediction_AAWA.mat');
        if(~exist(pred_file, 'file')); continue; end
        
        c = c+1;
        pred = load(pred_file);
        Nfolds = length(pred.acc);
        for f = 1:Nfolds
            data(c, b, f) = pred.pred_stats{f}(midx);
        end
    end
end
data = mean(data, 3);

%% plot
colormat = [200, 200, 200]./255;
f = figure('visible', 'off');
vio = violinplot(data, [], [], 'ViolinColor', colormat, 'ShowMean', true);
for i = 1:length(vio)
    vio(i).ViolinPlot.LineWidth = 2;
    vio(i).ScatterPlot.Marker = '.';
    vio(i).MedianPlot.SizeData = 18;
end
hold on
xlimit = get(gca, 'xlim');
plot(xlimit, [0 0], ':k');
hold off

pf = get(gcf, 'position');
set(gcf, 'position', [0 0 200+50*nbhvr 900])
set(gca, 'position', [0.45 0.4 0.6 0.5])

switch metric
case 'predictive_COD'
    y_label = 'Cross-validated predictive COD';
case 'corr'
    y_label = 'Cross-validated Pearson''s correlation';
otherwise
    error('Unknown metric')
end

yl = ylabel(y_label);
set(yl, 'fontsize', 16, 'linewidth', 2)
    
set(gca, 'xticklabel', ABCD_colloqs, 'fontsize', 16, 'linewidth', 2);
rotateXLabels( gca(), 45 );
set(gca, 'tickdir', 'out', 'box', 'off')

if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
outname = fullfile(outdir, [outstem ]);
export_fig(outname, '-png', '-nofontswap', '-a1');
set(gcf, 'color', 'w');
hgexport(f, outname)
close

rmpath(genpath( fullfile(getenv('HOME'), 'storage', 'from_HOME', 'code', 'plotting_functions')))
    
end