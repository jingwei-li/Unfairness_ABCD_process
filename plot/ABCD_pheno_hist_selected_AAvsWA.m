function ABCD_pheno_hist_selected_AAvsWA(csvfile, colloq_ls, AAWA_file, hist_dir, outstem)

addpath(genpath( fullfile(getenv('HOME'), 'storage', 'from_HOME', 'code', 'plotting_functions') ))

%% set default input arguments
ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');
if(~exist('csvfile', 'var') || isempty(csvfile))
    csvfile = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end
d = readtable(csvfile);
subj_hdr = 'subjectkey';

if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
[colloq_nm, ncolloq] = CBIG_text2cell(colloq_ls);

mkdir(hist_dir)
facecolor = [218,160,61; 97,98,71]./255;

%% load AAWA_file
% AAWA_file should contain variables: 
% (1) 'sel_AA', 'sel_WA', which are the selected subjects of each race group
% (2) 'bhvr_nm', behavioral names
load(AAWA_file)
assert(length(bhvr_nm) == ncolloq, 'Number of phenotypes is different from number of colloquial names.')

%% plot histograms
for b = 1:ncolloq
    % behavioral scores of each group
    AA = cat(2, selAA{b,:});
    WA = cat(1, selWA{b,:});
    [~, idxAA] = intersect(d.(subj_hdr), AA, 'stable');
    [~, idxWA] = intersect(d.(subj_hdr), WA, 'stable');
    mAA = d.(bhvr_nm{b})(idxAA);
    mWA = d.(bhvr_nm{b})(idxWA);
    
    % find bin edges
    alldata = [mAA; mWA];
    uq = unique(alldata);
    if(length(uq) > 50)
        nbin = 30;
    elseif(length(uq) > 20)
        nbin = 20;
    else
        nbin = length(uq);
    end
    h = histogram(alldata, nbin);
    E = h.BinEdges;
    binwidth = E(2) - E(1);
    close(gcf)
    
    % count frequencies per bin
    hc_WA = histcounts(mWA, E);
    hc_AA = histcounts(mAA, E);
    xloc = E(1:end-1) + diff(E)/2;
    
    % plot
    bh = bar(xloc, [hc_WA; hc_AA]');
    bh(1).FaceColor = facecolor(1,:);   bh(1).EdgeColor = 'none';
    bh(2).FaceColor = facecolor(2,:);   bh(2).EdgeColor = 'none';
    box off
    set(gcf, 'Position', [0 0 1600 600])
    
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
    if(binwidth<0.5 && binwidth>=0.01)
        xticklabels(sprintfc('%.2f', xloc))
    elseif(binwidth>10 || binwidth<0.01)
        %fprintf('%s: Keep original xticklabels\n', colloq_nm{b})
        xticklabels(sprintfc('%.2e', xloc))
        rotateXLabels( gca(), 30 );
    else
        xticklabels(sprintfc('%d', round(xloc)))
    end
    
    legend({'WA', 'AA'}, 'FontSize', 13, 'location', 'best')
    legend boxoff
    
    xloc_txt = xloc;
    text(xloc-binwidth/2 * 0.8, hc_WA + 0.05 * max([hc_AA hc_WA]), string(hc_WA), 'FontSize', 13)
    text(xloc, hc_AA + 0.05 * max([hc_AA hc_WA]), string(hc_AA), 'FontSize', 13)
    
    fname = fullfile(hist_dir, [colloq_nm{b} outstem '_WAvsAA']);
    [imageData, alpha] = export_fig([fname '.png'], '-png', '-nofontswap', '-a1');
    hgexport(gcf, fname)
    close(gcf)
end

rmpath(genpath( fullfile(getenv('HOME'), 'storage', 'from_HOME', 'code', 'plotting_functions') ))


end

