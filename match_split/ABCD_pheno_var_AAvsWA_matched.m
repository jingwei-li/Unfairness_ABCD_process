function ABCD_pheno_var_AAvsWA_matched(AAWA_file, full_subj_ls, pheno_csv, outname)

% ABCD_pheno_var_AAvsWA_matched(AAWA_file, full_subj_ls, pheno_csv)
%
% Conduct Levene's test to compare the variance of behavioral scores between 
% matched AA and WA.
%
% Inputs:
%   - AAWA_file
%     The .mat file containing selected matched AA-WA pairs for each site.
%     It is the output file of 'ABCD_match_and_split.m'.
%   - full_subj_ls
%     The list of all subject IDs (full path).
%   - pheno_csv
%     The csv file containing all the needed behavioral measures and demographics.
%     It is the output file of 'ABCD_read_all_measures.m'.
%   - outname
%     Full path of the output .mat file.
%
% Author: Jingwei Li

%% read input data
[full_subj, nsub] = CBIG_text2cell(full_subj_ls);
subj_hdr = 'subjectkey';
d = readtable(pheno_csv);
alpha = 0.05;

%% load AAWA_file
% AAWA_file should contain variables: 
% (1) 'sel_AA', 'sel_WA', which are the selected subjects of each race group
% (2) 'bhvr_nm', behavioral names
load(AAWA_file)
nbhvr = size(selAA, 1);

%% compare variances for each behavioral measure
p = zeros(nbhvr, 1);
stats = zeros(nbhvr, 1);
for b = 1:nbhvr
    % concatenate subject IDs across sites
    AA = cat(2, selAA{b,:});
    WA = cat(1, selWA{b,:});

    % extract corresponding behavioral scores
    [~, idxAA] = intersect(d.(subj_hdr), AA, 'stable');
    [~, idxWA] = intersect(d.(subj_hdr), WA, 'stable');
    mAA = d.(bhvr_nm{b})(idxAA);
    mWA = d.(bhvr_nm{b})(idxWA);

    % Levene's test
    [p(b), curr_stats] = vartestn([mWA; mAA], ...
        [ones(length(mWA),1); 2.*ones(length(mAA),1)], 'TestType', 'LeveneAbsolute');
    stats(b) = curr_stats.fstat;
end

%% FDR correction
H_idx = FDR(p(:), alpha);
H_FDR = zeros(size(p));
H_FDR(H_idx) = 1;
save(outname, 'p', 'stats', 'alpha', 'H_FDR')
    
end