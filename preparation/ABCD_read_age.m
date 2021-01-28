function [age, age_hdr] = ABCD_read_age(subj_list, race, dohist, hist_fname)

addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

age_csv = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/raw/documents/release2.0/ABCDstudyNDA/abcd_lt01.txt';
age_hdr = 'interview_age';
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

d = readtable(age_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');
age = cell(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        tmp_idx = tmp_idx & base_event;
        age(s) = d.(age_hdr)(tmp_idx);
    end
end

if(dohist==1)
    %% plot pure age distribution, without considering races
    age_plot = cellfun(@str2num, age);
    E = (min(age_plot)-1):2:(max(age_plot)+1);
    
    hc = histcounts(age_plot, E);
    bar(E(1:end-1) + diff(E)/2, hc');
    
    box off
    set(gcf, 'Position', [0 0 800 600])
    xloc = E(1:end-1)+diff(E)/2;
    xtl = sprintfc('%d', xloc);
    xloc_txt = E(1:end-1);
    text(xloc_txt, hc+50, string(hc), 'FontSize', 13)
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 13, 'TickDir','out')
    xticklabels(xtl)
    
    [outdir] = fileparts(hist_fname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close(gcf)
    
    %% plot AA/WA histograms across ages
    WA_filter = strcmp(race, '1');
    AA_filter = strcmp(race, '2');
    WA_age = age_plot(WA_filter);
    AA_age = age_plot(AA_filter);
    
    hc_WA = histcounts(WA_age, E);
    hc_AA = histcounts(AA_age, E);
    bar(E(1:end-1) + diff(E)/2, [hc_WA; hc_AA]')
    box off
    set(gcf, 'Position', [0 0 1000 600])
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 13, 'TickDir', 'out')
    legend({'WA', 'AA'}, 'FontSize', 13, 'location', 'best')
    legend boxoff
    
    text(xloc-0.4, hc_WA+15, string(hc_WA), 'FontSize', 13)
    text(xloc, hc_AA+15, string(hc_AA), 'FontSize', 13)
    xticklabels(xtl)
    
    [~, outbase, outext] = fileparts(hist_fname);
    fname2 = fullfile(outdir, [outbase '_WAvsAA' outext]);
    [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
    close(gcf)
end


rmpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))


end

