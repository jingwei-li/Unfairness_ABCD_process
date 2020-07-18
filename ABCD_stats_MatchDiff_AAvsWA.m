function [p_tt, H_tt_FDR, FDR_th] = ABCD_stats_MatchDiff_AAvsWA(sel_mAA, sel_mWA)

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

