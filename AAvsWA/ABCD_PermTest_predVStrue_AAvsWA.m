function ABCD_PermTest_predVStrue_AAvsWA(group_diff, bhvr_ls, outmat)

% ABCD_PermTest_predVStrue_AAvsWA(group_diff, bhvr_ls, outmat)
%
% Input
%   - group_diff
%     A .mat file contains the accuracy of each race group, and the
%     original & predicted scores of each group (full path).
%
%   - bhvr_ls (optional)
%     Behavior list (full path, text file). Use the behaviors which passed
%     predictability criteria.
%     Default: '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/behavior_list.txt'
%
%   - outmat
%     Output filename (full path).
%
% Author: Jingwei Li

alpha = 0.5;
nperm = 1000;
grpdif = load(group_diff);

if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', ...
        'scripts', 'lists', 'behavior_list.txt');
end
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);

%% permutation
rng('default'); rng(1, 'twister')
null_ydiff_AA = nan(nbhvr, size(grpdif.AA_test,2), nperm);
null_ydiff_WA = nan(nbhvr, size(grpdif.AA_test,2), nperm);
ydiff_AA = nan(nbhvr, size(grpdif.AA_test,2));  ydiff_WA = ydiff_AA;
for b = 1:nbhvr
    for f = 1:size(grpdif.AA_test, 2)
        ydiff_AA(b,f) = mean(grpdif.AA_pred{b,f} - grpdif.AA_test{b,f});
        ydiff_WA(b,f) = mean(grpdif.WA_pred{b,f} - grpdif.WA_test{b,f});

        for i = 1:nperm
            % indices of subjects to exchange ethnic/racial group labels
            idx = round(rand(length(grpdif.AA_test{b, f}), 1));
            % randomly shuffle the ordering of WA subjects
            shuffle = datasample(1:length(grpdif.WA_test{b,f}), length(grpdif.WA_test{b,f}), 'replace', false);
            
            null_yt_AA = grpdif.AA_test{b,f};
            null_yt_AA(idx==1) = grpdif.WA_test{b,f}(shuffle(idx==1));
            null_yp_AA = grpdif.AA_pred{b,f};
            null_yp_AA(idx==1) = grpdif.WA_pred{b,f}(shuffle(idx==1));
            
            null_yt_WA = grpdif.WA_test{b,f};
            null_yt_WA(shuffle(idx==1)) = grpdif.AA_test{b,f}(idx==1);
            null_yp_WA = grpdif.WA_pred{b,f};
            null_yp_WA(shuffle(idx==1)) = grpdif.AA_pred{b,f}(idx==1);
            
            null_ydiff_AA(b,f,i) = mean(null_yp_AA - null_yt_AA);
            null_ydiff_WA(b,f,i) = mean(null_yp_WA - null_yt_WA);
        end
    end
end

null_ydiff_diff = null_ydiff_WA - null_ydiff_AA;
null_ydiff_diff = squeeze(mean(null_ydiff_diff, 2));
avg_ydiff_diff = mean(ydiff_WA - ydiff_AA, 2);

p_perm = nan(nbhvr,1);
for b = 1:nbhvr
    p_perm(b) = length(find( null_ydiff_diff(b,:) - mean(null_ydiff_diff(b,:)) > abs(avg_ydiff_diff(b)) | ...
        null_ydiff_diff(b,:) - mean(null_ydiff_diff(b,:)) < -abs(avg_ydiff_diff(b)) )) / nperm;
end
p_perm(nbhvr+1) = length(find( mean(null_ydiff_diff,1) - mean(mean(null_ydiff_diff,1),2) > abs(mean(avg_ydiff_diff)) | ...
    mean(null_ydiff_diff,1) - mean(mean(null_ydiff_diff,1),2) < -abs(mean(avg_ydiff_diff)) )) / nperm;

H_perm_all = FDR(p_perm(:), alpha);
H_perm = setdiff(H_perm_all, nbhvr+1);
sig_diff_idx = sort(H_perm);
sig_diff_bhvr = bhvr_nm(sig_diff_idx);
    
save(outmat, 'p_perm', 'H_perm_all', 'H_perm', 'sig_diff_idx', 'sig_diff_bhvr');

    
end