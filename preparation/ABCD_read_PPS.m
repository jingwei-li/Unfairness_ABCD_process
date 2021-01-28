function [PPS, PPS_hdr, PPS_colloquial] = ABCD_read_PPS(subj_list, race, dohist, hist_dir, hist_fstem)

% [PPS, PPS_hdr, PPS_colloquial] = ABCD_read_PPS(subj_list, race, dohist, hist_dir, hist_fstem)
%
% Read and plot histogram of necessary mesures from Pediatric Psychosis Questionnaire 
%
% Example:
% [PPS, PPS_hdr, PPS_colloquial] = ABCD_read_PPS([], race, [], '~/storage/MyProject/fairAI/ABCD_race/figures/demo_hist', '_pass_rs');
% where "race" is obtained from 
% race = ABCD_read_race([], [], '~/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');

addpath(genpath( '~/storage/from_HOME/code/plotting_functions/'))

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

PPS_csv = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/raw/documents/release2.0/ABCDstudyNDA/abcd_mhy02.txt';
PPS_hdr = {'pps_y_ss_number', 'pps_y_ss_severity_score'};
PPS_colloquial = {'Total prodromal psychosis symptoms', 'Prodromal psychosis severity'};
for c = 1:length(PPS_colloquial)
    PPS_col_plot{c} = regexprep(PPS_colloquial{c}, ' +', '_');
end
subj_hdr = 'subjectkey';
event_hdr = 'eventname';

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

d = readtable(PPS_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');

% choose columns of selected PPS measures
PPS_read = [];
for c = 1:length(PPS_hdr)
    curr_PPS = d.(PPS_hdr{c});
    PPS_read = [PPS_read curr_PPS];
end

% select only the rows corresponding to required subjects
PPS = cell(nsub, length(PPS_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        tmp_idx = tmp_idx & base_event;
        PPS(s,:) = PPS_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, PPS);
PPS(empty_idx) = {'NaN'};
PPS = cellfun(@str2num, PPS);

if(dohist==1)
    binwidth = [1 5];
    nan_replace = [-2 -10];
    
    for c = 1:length(PPS_hdr)
        PPS_plot = PPS(:,c);
        % assign NaN to 10 (an invalid number for this task), so that #subjects without scores can be plotted
        PPS_plot(isnan(PPS_plot)) = nan_replace(c); 
        
        %% histogram across all subjects
        h = histogram(PPS_plot, 'binwidth', binwidth(c));
        box off
        set(gcf, 'Position', [0 0 1500 600])
        E = h.BinEdges;
        y = h.BinCounts; 
        start_idx = find(y~=0); start_idx = start_idx(2)-1;
        y(2:start_idx) = [];
        
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        xloc_txt = E(1:end-1); xloc_txt(2:start_idx) = [];
        text(xloc_txt, y+30, string(y), 'FontSize', 13)
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
        xticklabels(sprintfc('%d', [nan round(xloc(2:end))]))
        
        if(~exist(hist_dir, 'dir'))
            mkdir(hist_dir);
        end
        [imageData, alpha] = export_fig(fullfile(hist_dir, [PPS_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_PPS = PPS_plot(WA_filter);
        AA_PPS = PPS_plot(AA_filter);
        
        hc_WA = histcounts(WA_PPS, E); hc_WA(2:start_idx) = [];
        hc_AA = histcounts(AA_PPS, E); hc_AA(2:start_idx) = [];
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        bar(xloc, [hc_WA; hc_AA]')
        box off
        set(gcf, 'Position', [0 0 1500 600])
        
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
        xticklabels(sprintfc('%d', [nan round(xloc(2:end))]))
        
        legend({'WA', 'AA'}, 'FontSize', 13)
        legend boxoff
        
        xloc_txt = xloc([1 (start_idx+1):end]);
        text(xloc-binwidth(c)/2, hc_WA+30, string(hc_WA), 'FontSize', 13)
        text(xloc, hc_AA+30, string(hc_AA), 'FontSize', 13)
        
        fname2 = fullfile(hist_dir, [PPS_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end


rmpath(genpath( '~/storage/from_HOME/code/plotting_functions/'))


end

