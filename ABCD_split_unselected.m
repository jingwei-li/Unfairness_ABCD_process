function [fold_subj, fold_WA] = ABCD_split_unselected(subjects, site, sel_AA, sel_WA, fold_AA)

% [fold_subj] = ABCD_split_unselected(subjects, site, sel_AA, sel_WA, fold_AA)
%
% Inputs:
% - subjects: 1xN cell of strings, all subject IDs
% - site: Nx1 vector, site ID of all subjects
% - sel_AA: 1 x #site cell of strings, subject IDs of selected AA per site
% - sel_WA: 1 x #site cell of strings, subject IDs of selected WA per site
% - fold_AA: #fold x 1 cell of strings, subject IDs of selected AA in each fold
%

Nsites = length(sel_AA);
Nfolds = length(fold_AA);
uniq_sites = unique(site);

fold_subj = cell(Nfolds, 1);
fold_WA = cell(Nfolds, 1);
fold_sz = zeros(Nfolds, 1);
used_site = zeros(Nsites, 1);
for f = 1:Nfolds
    [~, AAidx] = intersect(subjects, fold_AA{f}');
    curr_sites = site(AAidx);
    curr_sites = unique(curr_sites);
    
    for st = 1:length(curr_sites)
        st_idx = uniq_sites == curr_sites(st);
        used_site(st_idx) = 1;
        
        % each fold contains selected AA and correspnding WA
        sel_AAWA = reshape(union(sel_AA{st_idx}, sel_WA{st_idx}), length(sel_AA{st_idx}) + length(sel_WA{st_idx}), 1);
        fold_subj{f} = [fold_subj{f}; sel_AAWA];
        fold_WA{f} = [fold_WA{f}; sel_WA{st_idx}];
        
        % remaining subjects belong to the same site are also in this fold
        subj_cr_st = subjects(site == curr_sites(st))';
        subj_cr_st = setdiff(subj_cr_st, sel_AAWA, 'stable');
        fold_subj{f} = [fold_subj{f}; subj_cr_st];
    end
    fold_sz(f) = length(fold_subj{f});
end

site_sz = zeros(Nsites, 1);
for st = 1:Nsites
    site_sz(st) = length(find(site == uniq_sites(st) ));
end
[~, sort_idx] = sort(site_sz, 'descend');

unused = find(used_site==0);
for st = 1:length(unused)
    curr_site = uniq_sites(sort_idx(st));
    
    [~, minf] = min(fold_sz);
    fold_subj{minf} = [fold_subj{minf}; subjects(site == curr_site)'];
    fold_sz(minf) = length(fold_subj{minf});
    
    used_site(unused(st)) = 1;
end


for f = 1:Nfolds
    fold_subj{f} = intersect(subjects', fold_subj{f}, 'stable');
end

end

