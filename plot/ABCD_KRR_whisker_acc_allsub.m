function ABCD_KRR_whisker_acc_allsub(model_dir, metric, outdir, outstem, bhvr_ls, colloq_ls)

% ABCD_KRR_whisker_acc_allsub(model_dir, metric, outdir, outname, bhvr_ls, colloq_ls)
%
% Plot the model performance in general (i.e. the accuracy on the whole test sets).

addpath(genpath( fullfile(getenv('HOME'), 'storage', 'from_HOME', 'code', 'plotting_functions')))
colormat = [200, 200, 200]./255;

proj_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race');
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(proj_dir, 'scripts', 'lists', 'behavior_list.txt');
end
if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(proj_dir, 'scripts', 'lists', 'colloquial_list.txt');
end

[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
colloq_nm = CBIG_text2cell(colloq_ls);
if(length(colloq_nm) ~= nbhvr)
    error('Length of behavioral names list and colloquial names list not equal.')
end

data = [];
[flag, msg] = system(['ls ' fullfile(model_dir, 'final_result_*.mat')]);
for b = 1:nbhvr
    if(flag==0)
        final_result = load(fullfile(model_dir, ['final_result_' bhvr_nm{b} '.mat']));
    else
        final_result = load(fullfile(model_dir, bhvr_nm{b}, ['final_result_' bhvr_nm{b} '.mat']));
    end
    data = [data final_result.optimal_stats.(metric)];
end
[~, idx] = sort(mean(data,1), 'descend');
data = data(:, idx);
bhvr_nm = bhvr_nm(idx);
colloq_nm = colloq_nm(idx);

%% plot
f = figure('visible', 'off');
aboxplot(data, 'colormap', colormat);
hold on
xlimit = get(gca, 'xlim');
plot(xlimit, [0 0], ':k');
hold off

pf = get(gcf, 'position');
set(gcf, 'position', [0 0 100+50*nbhvr 900])
set(gca, 'position', [0.35 0.4 0.6 0.5])

switch metric
    case 'corr'
        y_label = 'Cross-validated Pearson''s r';
    case 'predictive_COD'
        y_label = 'Cross-validated predictive COD';
    case 'COD'
        y_label = 'Cross-validated COD';
    otherwise
        error('Unknown metric')
end
yl = ylabel(y_label);
set(yl, 'fontsize', 16, 'linewidth', 2)

set(gca, 'xticklabel', colloq_nm, 'fontsize', 16, 'linewidth', 2);
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