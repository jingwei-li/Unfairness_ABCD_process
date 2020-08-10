function ABCD_KRR_scatter_AAvsWA(colloq_ls, pCOD_diff, outdir)

addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

%% figure configuration parameters setup
colorAA = [114 147 203] ./ 255;
colorWA = [132 186 91] ./ 255;
sz = 25;

%% parse input arguments
ls_dir = '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists';
if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
[colloq_nm, nbhvr] = CBIG_text2cell(colloq_ls);
% replace '_' with ' '
for b = 1:nbhvr
    curr_colloq = strsplit(colloq_nm{b}, '_');
    colloq_nm{b} = strjoin(curr_colloq, ' ');
end

grpdif = load(pCOD_diff);
CV = size(grpdif.AA_test, 2);

if(size(grpdif.AA_test, 1) ~= nbhvr)
    error('Numbers of behaviors in colloq_ls and group_diff are different')
end

%% scatter plot
corr_AA = nan(nbhvr, CV); corr_WA = corr_AA;
shift_sq_AA = nan(nbhvr, CV); shift_sq_WA = shift_sq_AA;
for b = 1:nbhvr
    if(~exist(fullfile(outdir, colloq_nm{b}), 'dir'))
        mkdir(fullfile(outdir, colloq_nm{b}))
    end
    fprintf('#%d behavior: %s\n', b, colloq_nm{b})

    for cv = 1:CV
        outname = fullfile(outdir, colloq_nm{b}, ['CV_' num2str(cv)]);
        if(exist([outname '.png'], 'file'))
            continue
        end
        
        f = figure('visible', 'off');
        scatter(grpdif.AA_test{b,cv}, grpdif.AA_pred{b,cv}, sz, colorAA, 'filled')
        hold on
        scatter(grpdif.WA_test{b,cv}, grpdif.WA_pred{b,cv}, sz, colorWA, 'filled')
        xli = get(gca, 'xlim');
        xpoints = xli(1):((xli(2)-xli(1))/5):xli(2);

        p = polyfit(grpdif.AA_test{b,cv}, grpdif.AA_pred{b,cv}, 1);
        r = polyval(p, xpoints);
        plot(xpoints, r, 'color', colorAA, 'LineWidth', 2)

        p = polyfit(grpdif.WA_test{b,cv}, grpdif.WA_pred{b,cv}, 1);
        r = polyval(p, xpoints);
        plot(xpoints, r, 'color', colorWA, 'LineWidth', 2)

        corr_AA(b,cv) = CBIG_corr(grpdif.AA_test{b,cv}, grpdif.AA_pred{b,cv});
        corr_WA(b,cv) = CBIG_corr(grpdif.WA_test{b,cv}, grpdif.WA_pred{b,cv});

        shift_sq_AA(b,cv) = (mean(grpdif.AA_test{b,cv} - grpdif.AA_pred{b,cv}))^2;
        shift_sq_WA(b,cv) = (mean(grpdif.WA_test{b,cv} - grpdif.WA_pred{b,cv}))^2;

        title_cell = {colloq_nm{b}, ...
            sprintf('corr-AA: %.2e, corr-WA: %.2e', corr_AA(b,cv), corr_WA(b,cv)), ...
            sprintf('pCOD-AA: %.2e, pCOD-WA: %.2e', grpdif.pCOD_AA(b,cv), grpdif.pCOD_WA(b,cv)), ...
            sprintf('mean-yt-AA: %.2e, mean-yt-WA: %.2e', ...
            mean(grpdif.AA_test{b,cv}), mean(grpdif.WA_test{b,cv})), ...
            sprintf('var-yt-AA: %.2e, mean-yt-WA: %.2e', ...
            var(grpdif.AA_test{b,cv},0), var(grpdif.WA_test{b,cv},0)), ...
            sprintf('(mean-yt-yp-diff)^2: AA: %.2e, WA: %.2e', shift_sq_AA(b,cv), shift_sq_WA(b,cv))};
        title(title_cell)
        set(gcf, 'position', [0 0 550 630])
        set(gca, 'linewidth', 2, 'fontsize', 13)
        xlabel('True scores', 'fontsize', 14)
        ylabel('Predicted scores', 'fontsize', 14)
        l = legend('AA', 'WA');
        set(l, 'fontsize', 13, 'Location', 'Best', 'box', 'off')

        export_fig(outname, '-png', '-nofontswap', '-a1');
        set(gcf, 'color', 'w')
        hgexport(f, outname)
        close
    end
end

rmpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

end

