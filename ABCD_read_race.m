function [race, race_hdr] = ABCD_read_race(subj_list, dohist, hist_fname)

% [race, race_hdr] = ABCD_read_race(subj_list, dohist, hist_fname)
%
% Example:
% race = ABCD_read_race([], [], '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');

addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

race_csv = '/mnt/eql/yeo12/data/ABCD/documents/release2.0/ABCDstudyNDA/acspsw03.txt';
race_hdr = 'race_ethnicity';
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

race_labels = {'', '1', '2', '3', '4', '5'};
race_grps = {'Empty', 'White', 'Black', 'Hispanic', 'Asian', 'Other'};
d = readtable(race_csv);
race = cell(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        race(s) = d.(race_hdr)(tmp_idx);
    end
end

if(dohist)
    uniq_race = unique(race);
    [~,idx] = intersect(race_labels, uniq_race);
    xtl = race_grps(idx);
    race_plot = race;
    for l = 1:length(race_labels)
        idx = strcmp(race, race_labels{l});
        race_plot(idx) = {l-1};
    end
    race_plot = cell2mat(race_plot);
    
    h = histogram(race_plot, 'BinWidth', 0.9999);
    box off
    set(gcf, 'Position', [0 0 600 600])
    E = h.BinEdges;
    y = h.BinCounts;
    xloc = E(1:end-1)+diff(E)/2;
    xloc_txt = E(1:end-1)+diff(E)/3;
    text(xloc_txt, y+100, string(y), 'FontSize', 13)
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 13, 'TickDir','out')
    xticklabels(xtl);
    
    [outdir] = fileparts(hist_fname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close
end

rmpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))


end

