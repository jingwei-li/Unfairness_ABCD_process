function [RAVLT, RAVLT_hdr, RAVLT_colloquial] = ABCD_read_RAVLT(subj_list, race, dohist, hist_dir, hist_fstem)

% [RAVLT, RAVLT_hdr, RAVLT_colloquial] = ABCD_read_RAVLT(subj_list, race, dohist, hist_dir, hist_fstem)
%
% Read necessary Rey Auditory Verbal Learning Test scores
%
% Example:
% RAVLT = ABCD_read_RAVLT([], race, [], '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist', '_pass_rs');
% where "race" is obtained from
% race = ABCD_read_race([], [], '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');

addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

RAVLT_csv = '/mnt/eql/yeo12/data/ABCD/documents/release2.0/ABCDstudyNDA/abcd_ps01.txt';
RAVLT_hdr = {'pea_ravlt_sd_trial_vi_tc', 'pea_ravlt_ld_trial_vii_tc'};
RAVLT_colloquial = {'Short delay recall', 'Long delay recall'};
for c = 1:length(RAVLT_colloquial)
    RAVLT_col_plot{c} = regexprep(RAVLT_colloquial{c}, ' +', '_');
end
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

d = readtable(RAVLT_csv);

% choose columns of selected RAVLT measures
RAVLT_read = [];
for c = 1:length(RAVLT_hdr)
    curr_RAVLT = d.(RAVLT_hdr{c});
    RAVLT_read = [RAVLT_read curr_RAVLT];
end

% select only the rows corresponding to required subjects
RAVLT = cell(nsub, length(RAVLT_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        RAVLT(s,:) = RAVLT_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, RAVLT);
RAVLT(empty_idx) = {'NaN'};
RAVLT = cellfun(@str2num, RAVLT);

if(dohist==1)
    for c = 1:length(RAVLT_hdr)
        RAVLT_plot = RAVLT(:,c);
        % assign NaN to -5 (an invalid number for this task), so that #subjects without scores can be plotted
        RAVLT_plot(isnan(RAVLT_plot)) = -5; 
        
        %% histogram across all subjects
        h = histogram(RAVLT_plot);
        box off
        set(gcf, 'Position', [0 0 800 600])
        E = h.BinEdges;
        y = h.BinCounts; y(2:5) = [];
        xloc = unique(RAVLT_plot);
        xloc_txt = E(1:end-1); xloc_txt(2:5) = [];
        text(xloc_txt, y+20, string(y), 'FontSize', 13)
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
        xticklabels(sprintfc('%d', [nan; xloc(2:end)]))
        
        if(~exist(hist_dir, 'dir'))
            mkdir(hist_dir);
        end
        [imageData, alpha] = export_fig(fullfile(hist_dir, [RAVLT_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_RAVLT = RAVLT_plot(WA_filter);
        AA_RAVLT = RAVLT_plot(AA_filter);
        
        hc_WA = histcounts(WA_RAVLT, E); hc_WA(2:5) = [];
        hc_AA = histcounts(AA_RAVLT, E); hc_AA(2:5) = [];
        xloc = unique([WA_RAVLT; AA_RAVLT]);
        bar(xloc, [hc_WA; hc_AA]')
        box off
        set(gcf, 'Position', [0 0 1300 600])
        
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
        xticklabels(sprintfc('%d', [nan; xloc(2:end)]))
        
        legend({'WA', 'AA'}, 'FontSize', 13)
        legend boxoff
        
        xloc_txt = xloc([1 6:end]);
        text(xloc-0.4, hc_WA+10, string(hc_WA), 'FontSize', 13)
        text(xloc, hc_AA+10, string(hc_AA), 'FontSize', 13)
        
        fname2 = fullfile(hist_dir, [RAVLT_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end

rmpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

end

