function income = ABCD_read_income(subj_list, race, dohist, hist_fname)

% Example:
% income = ABCD_read_income([], race, [], '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/income_pass_rs.png');
% where "race" is obtained from
% race = ABCD_read_race([], [], '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');

addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

income_csv = '/mnt/eql/yeo12/data/ABCD/documents/release2.0/ABCDstudyNDA/pdem02.txt';
income_hdr = 'demo_comb_income_v2';
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

d = readtable(income_csv);
income = cell(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        income(s) = d.(income_hdr)(tmp_idx);
    end
end
empty_idx = cellfun(@isempty, income);
income{empty_idx} = '999';

if(dohist==1)
    % 1= Less than $5,000; 2=$5,000 through $11,999; 3=$12,000 through $15,999; 4=$16,000 through $24,999; 
    % 5=$25,000 through $34,999; 6=$35,000 through $49,999; 7=$50,000 through $74,999; 
    % 8= $75,000 through $99,999; 9=$100,000 through $199,999; 10=$200,000 and greater. 
    % 999 = Don't know; 777 = Refuse to answer
    xtl = {'Refuse ans', 'Unknown', '<$5k', '$5k-12k', '$12k-16k', '$16k-25k', ...
        '$25k-35k', '$35k-50k', '$50k-75k', '$75k-100k', '$100k-200k', '>=$200k'};
    
    %% plot household income of all subjects
    income_plot = cellfun(@str2num, income);
    income_plot(income_plot==777) = -0.99;
    income_plot(income_plot==999) = 0;
    
    h = histogram(income_plot, 'BinWidth', 0.9999);
    box off
    set(gcf, 'Position', [0 0 1400 600])
    E = h.BinEdges;
    y = h.BinCounts;
    xloc = E(1:end-1)+diff(E)/2;
    xloc_txt = E(1:end-1) + diff(E)/3;
    text(xloc_txt, y+50, string(y), 'FontSize', 13)
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
    xticklabels(xtl)
    
    [outdir] = fileparts(hist_fname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close(gcf)
    
    %% plot household income for WA/AA separately
    WA_filter = strcmp(race, '1');
    AA_filter = strcmp(race, '2');
    WA_income = income_plot(WA_filter);
    AA_income = income_plot(AA_filter);
    
    hc_WA = histcounts(WA_income, E);
    hc_AA = histcounts(AA_income, E);
    bar(E(1:end-1), [hc_WA; hc_AA]')
    box off
    set(gcf, 'Position', [0 0 1300 600])
    xloc_txt = min(round(income_plot)):1:max(income_plot);
    set(gca, 'xtick', xloc_txt, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
    xticklabels(xtl)
    
    legend({'WA', 'AA'}, 'FontSize', 13)
    legend boxoff
    
    text(xloc_txt-0.4, hc_WA+30, string(hc_WA), 'FontSize', 13)
    text(xloc_txt, hc_AA+30, string(hc_AA), 'FontSize', 13)
    
    [~, outbase, outext] = fileparts(hist_fname);
    fname2 = fullfile(outdir, [outbase '_WAvsAA' outext]);
    [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
    close(gcf)
end

rmpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

end

