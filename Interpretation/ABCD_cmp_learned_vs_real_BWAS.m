function [simi, dst] = ABCD_cmp_learned_vs_real_BWAS(learned_BWAS, real_BWAS, outmat, fig_dir)

% ABCD_cmp_learned_vs_real_BWAS(learned_BWAS, real_BWAS, outmat, fig_dir)
%
% 

load(learned_BWAS, 'avg_learned_cov')
load(real_BWAS, 'avg_cov_testAA', 'avg_cov_testWA')
nrois = size(avg_learned_cov, 1);
tril_mtx = ones(nrois, nrois);
tril_mtx = tril(tril_mtx, -1);

nbhvr = size(avg_learned_cov, 3);
avg_learned_cov = reshape(avg_learned_cov, nrois^2, nbhvr);
avg_cov_testAA = reshape(avg_cov_testAA, nrois^2, nbhvr);
avg_cov_testWA = reshape(avg_cov_testWA, nrois^2, nbhvr);

avg_learned_cov = avg_learned_cov(tril_mtx(:)==1, :);
avg_cov_testAA = avg_cov_testAA(tril_mtx(:)==1, :);
avg_cov_testWA = avg_cov_testWA(tril_mtx(:)==1, :);

sim_AA = zeros(nbhvr, 1);   sim_WA = zeros(nbhvr, 1);
dist_AA = zeros(nbhvr, 1);  dist_WA = zeros(nbhvr, 1);
for b = 1:nbhvr
    sim_AA(b) = CBIG_corr(avg_learned_cov(:,b), avg_cov_testAA(:,b));
    sim_WA(b) = CBIG_corr(avg_learned_cov(:,b), avg_cov_testWA(:,b));

    dist_AA(b) = sqrt(sum((avg_learned_cov(:,b) - avg_cov_testAA(:,b)).^2));
    dist_WA(b) = sqrt(sum((avg_learned_cov(:,b) - avg_cov_testWA(:,b)).^2));
end
simi = [sim_AA sim_WA];
dst = [dist_AA dist_WA];

    
end