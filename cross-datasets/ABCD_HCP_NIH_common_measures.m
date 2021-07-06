function ABCD_HCP_NIH_common_measures(fig_dir, ABCD_subj_ls, HCP_subj_ls, ABCD_csv, HCP_unrstr_csv)

% ABCD_HCP_NIH_common_measures()
%
% For the common NIH measures used in both datasets, check if the scores are
% on the same scale

ABCD_bhvr = {
    'nihtbx_flanker_uncorrected', ...
    'nihtbx_list_uncorrected', ...
    'nihtbx_cardsort_uncorrected', ...
    'nihtbx_reading_uncorrected', ...
    'nihtbx_pattern_uncorrected', ...
    'nihtbx_picture_uncorrected', ...
    'nihtbx_picvocab_uncorrected', ...
};
HCP_bhvr = {
    'Flanker_Unadj', ...
    'ListSort_Unadj', ...
    'CardSort_Unadj', ...
    'ReadEng_Unadj', ...
    'ProcSpeed_Unadj', ...
    'PicSeq_Unadj', ...
    'PicVocab_Unadj', ...
};

proj_dir = '/home/jingweil/storage/MyProject/fairAI/';
if(~exist('ABCD_subj_ls', 'var') || isempty(ABCD_subj_ls))
    ABCD_subj_ls = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'subjects_pass_rs_pass_pheno.txt');
end
if(~exist('HCP_subj_ls', 'var') || isempty(HCP_subj_ls))
    HCP_subj_ls = fullfile(proj_dir, 'HCP_race', 'scripts', 'lists', 'subjects_wIncome_948.txt');
end
if(~exist('ABCD_csv', 'var') || isempty(ABCD_csv))
    ABCD_csv = fullfile(proj_dir, 'ABCD_race', 'scripts', 'lists', 'phenotypes_pass_rs.txt');
end
if(~exist('HCP_unrstr_csv', 'var') || isempty(HCP_unrstr_csv))
    HCP_unrstr_csv = '/mnt/isilon/CSC1/Yeolab/Data/HCP/S1200/scripts/subject_measures/unrestricted_jingweili_12_7_2017_21_0_16_NEO_A_corrected.csv';
end

[ABCD_subj, N_abcd] = CBIG_text2cell(ABCD_subj_ls);
[HCP_subj, N_hcp] = CBIG_text2cell(HCP_subj_ls);
abcd_subj_hdr = 'subjectkey';
hcp_subj_hdr = 'Subject';
d_abcd = readtable(ABCD_csv);
d_hcp = readtable(HCP_unrstr_csv);

for b = 1:length(ABCD_bhvr)
    curr_ABCD = d_abcd.(ABCD_bhvr{b});
    [~, ~, abcd_idx] = intersect(ABCD_subj, d_abcd.(abcd_subj_hdr), 'stable');
    curr_ABCD = curr_ABCD(abcd_idx);

    curr_HCP = d_hcp.(HCP_bhvr{b});
    [~, ~, hcp_idx] = intersect(HCP_subj, cellfun(@num2str, ...
        mat2cell(d_hcp.(hcp_subj_hdr), ones(length(d_hcp.(hcp_subj_hdr)), 1), 1), 'UniformOutput', false), 'stable');
    curr_HCP = curr_HCP(hcp_idx);

    perc_ABCD = length(find(curr_ABCD > min(curr_HCP) & curr_ABCD < max(curr_HCP))) / N_abcd;
    perc_HCP = length(find(curr_HCP > min(curr_ABCD) & curr_HCP < max(curr_ABCD))) / N_hcp;

    fprintf('[ABCD]: %s, min: %f, max: %f\n', ABCD_bhvr{b}, min(curr_ABCD), max(curr_ABCD));
    fprintf('[HCP]: %s, min: %f, max: %f\n', HCP_bhvr{b}, min(curr_HCP), max(curr_HCP));
    fprintf('Percentage of overlap in ABCD: %f\n', perc_ABCD)
    fprintf('Percentage of overlap in HCP: %f\n', perc_HCP)
    fprintf('--------------------------\n')

    hist_fname = fullfile(fig_dir, [ABCD_bhvr{b} '_' HCP_bhvr{b} '.png']);
    %hist_2grp(curr_ABCD, curr_HCP, ABCD_bhvr{b}, HCP_bhvr{b}, hist_fname);
end
    
end


function hist_2grp(ABCD_scores, HCP_scores, ABCD_bhvr, HCP_bhvr, hist_fname)
    
    lb = min(min(ABCD_scores), min(HCP_scores));
    ub = max(max(ABCD_scores), max(ABCD_scores));
    n_uniq = length(unique([ABCD_scores; HCP_scores]));
    if(n_uniq<50)
        nbins = n_uniq;
    else
        nbins = 50;
    end
    binwidth = round((ub - lb) / nbins);
    E = (lb - binwidth/2):binwidth:(ub + binwidth/2);
    hc_abcd = histcounts(ABCD_scores, E);
    hc_hcp = histcounts(HCP_scores, E);
    xloc = E(1:end-1)+diff(E)/2;
    bar(xloc, [hc_abcd; hc_hcp]');

    box off
    set(gcf, 'Position', [0 0 1000 600])
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 11, 'TickDir', 'out')
    legend({['ABCD: ' strjoin(strsplit(ABCD_bhvr, '_'), ' ')], ...
        ['HCP: ' strjoin(strsplit(HCP_bhvr, '_'), ' ')]}, 'FontSize', 11, 'location', 'best')
    legend boxoff

    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close(gcf)
end