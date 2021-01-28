function ICV = ABCD_read_ICV(subj_list, race, dohist, hist_fname)

addpath(genpath( '~/storage/from_HOME/code/plotting_functions/'))

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

FSdir = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/process/y0/recon_all';

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt';
end
[subjects, nsub] = CBIG_text2cell(subj_list);

ICV = nan(nsub,1);
for s = 1:nsub
    FSstats = CBIG_text2cell(fullfile(FSdir, subjects{s}, 'stats', 'aseg.stats'));
    ICV_row = contains(FSstats, 'Intracranial Volume');
    if(exist('ICV_row', 'var') && any(ICV_row==1))
        ICV_row = FSstats{ICV_row==1};
        row_split = strsplit(ICV_row, ',');
        ICV(s) = str2double(row_split{4});
    end
end

if(dohist==1)
    binwidth = 5e4;
    E = (min(ICV)-binwidth/2):binwidth:(max(ICV)+binwidth/2);
    hc = histcounts(ICV, E);
    bar(E(1:end-1) + diff(E)/2, hc');
    
    box off
    set(gcf, 'Position', [0 0 800 600])
    xloc = E(1:end-1)+diff(E)/2; 
    xloc1 = xloc(1:4:length(xloc));
    xtl = sprintfc('%.2e', xloc1);
    xloc_txt = E(1:end-1);  
    text(xloc_txt, hc+50, string(hc), 'FontSize', 13)
    set(gca, 'xtick', xloc1, 'linewidth', 2, 'fontsize', 13, 'TickDir','out')
    xticklabels(xtl)
    
    [outdir] = fileparts(hist_fname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close(gcf)
    
    %% plot AA/WA histograms across peduc_avg
    WA_filter = strcmp(race, '1');
    AA_filter = strcmp(race, '2');
    WA_ICV = ICV(WA_filter);
    AA_ICV = ICV(AA_filter);
    
    hc_WA = histcounts(WA_ICV, E);
    hc_AA = histcounts(AA_ICV, E);
    bar(E(1:end-1) + diff(E)/2, [hc_WA; hc_AA]')
    box off
    set(gcf, 'Position', [0 0 1000 600])
    set(gca, 'xtick', xloc1, 'linewidth', 2, 'fontsize', 13, 'TickDir', 'out')
    legend({'WA', 'AA'}, 'FontSize', 13, 'location', 'best')
    legend boxoff
    
    text(xloc-binwidth/2, hc_WA+15, string(hc_WA), 'FontSize', 13)
    text(xloc, hc_AA+15, string(hc_AA), 'FontSize', 13)
    xticklabels(xtl)
    
    [~, outbase, outext] = fileparts(hist_fname);
    fname2 = fullfile(outdir, [outbase '_WAvsAA' outext]);
    [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
    close(gcf)
end

rmpath(genpath( '~/storage/from_HOME/code/plotting_functions/'))

end

