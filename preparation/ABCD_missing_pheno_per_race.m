function [nan_per_eth, nonnan_per_eth] = ABCD_missing_pheno_per_race(pheno_csv, bhvr_ls)

% [nan_per_eth, nonnan_per_eth] = ABCD_missing_pheno_per_race(pheno_csv, bhvr_ls)
%
% Check if the minorities have more missing values for behavioral measures, compared to 
% white Americans.
% 
% Inputs:
%   - pheno_csv
%     The csv file containing the needed behavioral measures for all subjects (who survived 
%     imaging quality control). It is collecte by 'ABCD_read_all_measures.m'.
%   - bhvr_ls
%     Full path of the list containing the behavioral measures of interest.
% 
% Outputs:
%   - nan_per_race
%     A #behavior x #ethnicities matrix. Entry (i,j) represents the number of missing values 
%     for the i-th behavioral measure in the j-th racial/ethnic group.
%   - nonnan_per_race
%     A #behavior x #ethnicities matrix. Entry (i,j) represents the number of subjects who 
%     have been recorded for the i-th behavioral measures in the j-th racial/ethnic group.
%
% Author: Jingwei Li

d = readtable(pheno_csv);
race_hdr = 'race';
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

%% collect the indicies of subjects belongning to each racial/ethnic group
uq_race = unique(d.(race_hdr));
uq_race(isnan(uq_race)) = [];
idx_per_race = cell(length(uq_race), 1);
for r = 1:length(uq_race)
    idx_per_race{r} = find(d.(race_hdr) == uq_race(r));
end

%% count subjects w/ or w/o behavioral scores per group
nan_per_eth = zeros(nbhvr, length(uq_race));
nonnan_per_eth = zeros(nbhvr, length(uq_race));
for b = 1:nbhvr
    pheno = d.(bhvr_nm{b});
    nan_idx = find(isnan(pheno));
    nonnan_idx = find(~isnan(pheno));

    for r = 1:length(uq_race)
        nan_per_eth(b,r) = length(intersect(idx_per_race{r}, nan_idx));
        nonnan_per_eth(b,r) = length(intersect(idx_per_race{r}, nonnan_idx));
    end
end
    
end