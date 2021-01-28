function [site, site_hdr] = ABCD_read_site(subj_list, race, dohist, hist_fname)

% [site, site_hdr] = ABCD_read_site(subj_list, race, dohist, hist_fname)
%
% Example:
% site = ABCD_read_site([], race, [], '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/site_pass_rs.png');
% "race" is calculated using 
% race = ABCD_read_race([], [], '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');

addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

site_csv = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/raw/documents/release2.0/ABCDstudyNDA/abcd_lt01.txt';
site_hdr = 'site_id_l';
event_hdr = 'eventname';
subj_hdr = 'subjectkey';

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt';
end
[subjects, nsub] = CBIG_text2cell(subj_list);

% format in subj_list: e.g. NDARINV007W6H7B
% format in csv: e.g. "NDAR_INV007W6H7B"
subjects_csv = cell(nsub, 1);
for s = 1:nsub
    subjects_csv{s} = [subjects{s}(1:4) '_' subjects{s}(5:end)];
end

d = readtable(site_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');
site = cell(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        tmp_idx = tmp_idx & base_event;
        site(s) = d.(site_hdr)(tmp_idx);
    end
end

if(dohist==1)
    %% plot pure site distribution, without considering raes
    site_plot = site;
    for i = 1:length(site_plot)
        site_plot{i} = str2num(site_plot{i}(5:end));
    end
    site_plot = cell2mat(site_plot);
    xtl = sort(unique(site_plot));
    xtl = sprintfc('%d', xtl);
    
    h = histogram(site_plot, 'BinWidth', 0.9999);
    box off
    set(gcf, 'Position', [0 0 800 600])
    E = h.BinEdges;
    y = h.BinCounts;
    xloc = E(1:end-1)+diff(E)/2;
    xloc_txt = E(1:end-1);
    text(xloc_txt, y+50, string(y), 'FontSize', 13)
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 13, 'TickDir','out')
    xticklabels(xtl)
    
    [outdir] = fileparts(hist_fname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close(gcf)
    
    %% plot AA/WA histograms across sites
    WA_filter = strcmp(race, '1');
    AA_filter = strcmp(race, '2');
    WA_site = site_plot(WA_filter);
    AA_site = site_plot(AA_filter);
    
    hc_WA = histcounts(WA_site, E);
    hc_AA = histcounts(AA_site, E);
    bar(E(1:end-1), [hc_WA; hc_AA]')
    box off
    set(gcf, 'Position', [0 0 1200 600])
    set(gca, 'linewidth', 2, 'fontsize', 13, 'TickDir', 'out')
    legend({'WA', 'AA'}, 'FontSize', 13)
    legend boxoff
    
    xloc_txt = min(site_plot):1:max(site_plot);
    text(xloc_txt-0.4, hc_WA+15, string(hc_WA), 'FontSize', 13)
    text(xloc_txt, hc_AA+15, string(hc_AA), 'FontSize', 13)
    
    [~, outbase, outext] = fileparts(hist_fname);
    fname2 = fullfile(outdir, [outbase '_WAvsAA' outext]);
    [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
    close(gcf)
end


end

