function ABCD_whisker_AAvsWA(bhvr_ls, colloq_ls, group_diff, perm_fname, metric, outdir, outstem)

% ABCD_whisker_AAvsWA(bhvr_ls, colloq_ls, group_diff, perm_fname, metric, outdir, outstem)
% 
% Input:
%   - bhvr_ls (optional)
%     Behavior list (full path, text file). Use the behaviors which passed
%     predictability criteria. Default: 
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
% 
%   - colloq_ls (optional)
%     List of behaviors' colloquial names (full path). Use the same
%     behaviors as in bhvr_ls. Default: 
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/colloquial_list.txt'
%
%   - group_diff
%     A .mat file contains the accuracy of each race group, and the
%     original & predicted scores of each group (full path).
%
%   - perm_fname
%     A .mat file contains the permutation testing results of AA-WA
%     accuracy difference (full path).
%
%   - metric
%     Type of accuracy metric. For now only support 'predictive_COD', 
%     'corr'.
%
%   - outdir
%     Output directory (full path).
%
%   - outstem
%     A string to be atttached onto the output filenames, e.g.
%     '_pass_rs_pass_pheno'.
%

addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

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
    otherwise
        error('Unknown metric.')
end

colormat = [114 147 203; 132 186 91; 211 94 96]./255;
legend1 = 'AA';
legend2 = 'Matched WA';
legend3 = 'Difference';

%% parse input arguments, collect data
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
colloq_nm = CBIG_text2cell(colloq_ls);
% replace '_' with ' '
for b = 1:nbhvr
    curr_colloq = strsplit(colloq_nm{b}, '_');
    colloq_nm{b} = strjoin(curr_colloq, ' ');
end

grpdif = load(group_diff);
load(perm_fname);

CV = size(grpdif.(AA_acc), 2);

acc_diff = grpdif.(WA_acc)' - grpdif.(AA_acc)';
alldata = cat(1, reshape(grpdif.(AA_acc)', 1, CV, nbhvr), reshape(grpdif.(WA_acc)', 1, CV, nbhvr), ...
    reshape(acc_diff, 1, CV, nbhvr));
[~, idx] = sort(mean(acc_diff, 1), 'descend');
data_sort = alldata(:,:,idx);
bhvr_nm_sort = bhvr_nm(idx);
colloq_nm_sort = colloq_nm(idx);

% get indices of significant behaviors
[~, IA, IB] = intersect(idx, sig_diff_idx, 'stable');

if(~exist(outdir, 'dir'))
    mkdir(outdir)
end

%% plot for each behavior
f = figure('visible', 'off');
aboxplot(data_sort, 'colormap', colormat)
hold on
xlimit = get(gca, 'xlim');
plot(xlimit, [0 0], ':k');
hold off

pf = get(gcf, 'position');
set(gcf, 'position', [0 0 100+90*nbhvr 900])
set(gca, 'position', [0.35 0.4 0.6 0.5])

ylm = get(gca, 'ylim');
if(ylm(1)<-1)
    ylm(1) = -1;
    warning('There are behaviors with accuracy lower than -1.')
end
set(gca, 'ylim', ylm, 'ytick', [ylm(1):0.2:ylm(2)]);
yl = ylabel(y_label);
set(yl, 'fontsize', 16, 'linewidth', 2)

l = legend(legend1, legend2, legend3);
set(l, 'fontsize', 12, 'linewidth', 2, 'location', 'best', 'box', 'off')

set(gca, 'xticklabel', colloq_nm_sort, 'fontsize', 16, 'linewidth', 2);
rotateXLabels( gca(), 45 );
set(gca, 'tickdir', 'out', 'box', 'off')

% plot * on the significant behaviors
if(~isempty(IA))
    ylimvals = get(gca, 'ylim');
    ylimvals_new = ylimvals;
    ylimvals_new(2) = (ylimvals(2)-ylimvals(1))*0.08 + ylimvals(2);
    set(gca, 'ylim', ylimvals_new)
    lp = get(l, 'position');
    if(lp(2)>0.5)
        lp(2) = lp(2) - 0.03;
    end
    set(l, 'position', lp)
    text(IA, repmat(ylimvals(2), size(IA)), '*');
end

outname = fullfile(outdir, [outstem ]);
export_fig(outname, '-png', '-nofontswap', '-a1');
set(gcf, 'color', 'w');
hgexport(f, outname)
close

%% plot the average
avg_data = mean(data_sort, 3)';
f = figure('visible', 'off');
aboxplot(avg_data, 'colormap', colormat, 'colorrev', 1);
hold on
xlimit = get(gca, 'xlim');
plot(xlimit, [0 0], ':k');
hold off

pf = get(gcf, 'position');
set(gcf, 'position', [0 0 300 800]);
set(gca, 'position', [0.3 0.3 0.6 0.6])
ylm = get(gca, 'ylim');
set(gca, 'ytick', ylm(1):0.2:ylm(2))
yl = ylabel(y_label_avg, 'fontsize', 16, 'linewidth', 2);

set(gca, 'XTickLabel', {legend1, legend2, legend3}, 'fontsize', 16, 'linewidth', 2)
rotateXLabels( gca(), 45 );
set(gca, 'tickdir', 'out', 'box', 'off');

outname = fullfile(outdir, ['Mean_' outstem]);
export_fig(outname, '-png', '-nofontswap', '-a1');
set(gcf, 'color', 'w');
hgexport(f, outname)
close

rmpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

end

