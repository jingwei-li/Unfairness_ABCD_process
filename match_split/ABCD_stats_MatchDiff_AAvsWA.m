function [p_tt, H_tt_FDR, FDR_th] = ABCD_stats_MatchDiff_AAvsWA(sel_mAA, sel_mWA)

% [p_tt, H_tt_FDR, FDR_th] = ABCD_stats_MatchDiff_AAvsWA(sel_mAA, sel_mWA)
%
% Two-sample t test to check whether the confounding variables and behavioral scores
% were still significantly different between matched AA and WA.
%
% Inputs:
% - sel_mAA
%   A cell array with a length: number of behavioral measures. Each entry is a 
%   matrix (dimension: #confounds+1 x #AA) containing the confounds and the scores of
%   a certain behavioral measure of all matched AA concatenated across sites.
%   It is the output of `ABCD_match_WAtoAA_within_site.m`.
%
% - sel_mWA
%   A cell array with a length: number of behavioral measures. Each entry is a 
%   matrix (dimension: #confounds+1 x #WA) containing the confounds and the scores of
%   a certain behavioral measure of all matched WA concatenated across sites.
%   It is the output of `ABCD_match_WAtoAA_within_site.m`.
%
% Outputs:
% - p_tt
%   A #behaviors x #matched_variables (i.e. #confounds+1) matrix. Each entry is the
%   p value resulted from t test, indicating the significancy of the difference
%   of a certain variable to be matched between the selected AA and AA, for a certain 
%   behavioral measure.
%
% - H_tt_FDR
%   A #behaviors x #matched_variables (i.e. #confounds+1) matrix with values of 0 or 1. 
%   Each entry indicates whether the null hypothesis (variable to be matched do not 
%   differ between AA and WA) is rejected or not (false discovery rate < `FDR_th`).
%
% - FDR_th
%   Threshold used for FDR correction. 
%
% Author: Jingwei Li

nbehav = length(sel_mAA);
nmetric = size(sel_mAA{1}, 2);
alpha = 0.05;

p_tt = zeros(nbehav, nmetric);
for b = 1:nbehav
    for m = 1:nmetric
        [~, p_tt(b,m)] = ttest2(sel_mAA{b}(:,m), sel_mWA{b}(:,m), 'Vartype', 'unequal');
    end
end

[H_tt_FDR_idx, FDR_th] = FDR(p_tt(:), alpha);
H_tt_FDR = zeros(size(p_tt));
H_tt_FDR(H_tt_FDR_idx) = 1;

end

