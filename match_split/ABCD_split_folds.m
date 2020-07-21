function [fold_subj, fold_sites] = ABCD_split_folds(subjects, site, fam_id, Nfolds)

if(~exist('Nfolds', 'var') || isempty(Nfolds))
    Nfolds = 10;
end

N = length(subjects);
N_fold = round(N / Nfolds);

if(size(subjects,2) > 1)
    subjects = subjects';
end

%% check if the subjects from the same family spread across sites
uniq_fam = unique(fam_id);
sites_same_fam = [];
fam_sz = zeros(length(uniq_fam), 1);
for f = 1:length(uniq_fam)
    idx = find(fam_id == uniq_fam(f));
    fam_sz(f) = length(idx);
    if(fam_sz(f)>1)
        fam_site = site(idx);
        if(length(unique(fam_site)) > 1)
            sites_same_fam = [site_same_fam; {unique(fam_site)}];
            str = strjoin(sprintfc('%d', unique(fam_site)), ',');
            warning('Sites %s must be put into the same fold.', str);
        end
    end
end

%% if there are subjects from the same family across different sites, combine those sites
if(~isempty(sites_same_fam))
    restrictions = cat(1, sites_same_fam{:});
    for i = 1:length(restrictoins)
        idx = find(cellfun(@(x)(ismember(restrictions(i), x)), sites_same_fam));
        sites_same_fam{min(idx)} = cat(1, sites_same_fam{idx});
        sites_same_fam(idx(2:end)) = [];
    end
    
    for i = 1:length(sites_same_fam)
        for l = 1:length(sites_same_fam{i})
            site(site == sites_same_fam{i}(l)) = min(sites_same_fam{i});
        end
    end
end

%% deal with too large individual sites
uniq_site = unique(site);
site_sz = zeros(length(uniq_site), 1);
site_subj = cell(length(uniq_site), 1);
for i = 1:length(site_sz)
    site_sz(i) = length(find(site == uniq_site(i) ));
    site_subj{i} = subjects(site == uniq_site(i));
end

exceed_sites = site_sz > N_fold;
N_exceeded = sum(site_sz(exceed_sites)) - N_fold * sum(double(exceed_sites));
if(N_exceeded >= N_fold)
    % if the large sites take up the size of an entire fold, adjust the
    % upper-bound size of average fold
    N_fold = round( (N - sum(site_sz(exceed_sites))) / (Nfolds - sum(double(exceed_sites))) );
end


%%
% [~, sort_idx] = sort(site_sz, 'descend');

rng('default')
rng(1000)
sort_idx = datasample(1:1:length(site_sz), length(site_sz), 'Replace', false);

used_site = zeros(length(uniq_site), 1);
fold_sz = zeros(Nfolds,1);
fold_sites = cell(Nfolds, 1);
for f = 1:Nfolds
    for s = 1:length(uniq_site)
        if(used_site(sort_idx(s))==0)
            if(fold_sz(f) + site_sz(sort_idx(s)) <= N_fold || (fold_sz(f)==0 && site_sz(sort_idx(s)) > N_fold) )
                fold_sz(f) = fold_sz(f) + site_sz(sort_idx(s));
                fold_sites{f} = [fold_sites{f}; uniq_site(sort_idx(s))];
                used_site(sort_idx(s)) = 1;
            end
        end
    end
end

sz_dif = N_fold - fold_sz;
for s = 1:length(uniq_site)
    if(used_site(sort_idx(s)) == 0)
        [~,max_f] = max(sz_dif);
        fold_sz(max_f) = fold_sz(max_f) + site_sz(sort_idx(s));
        fold_sites{max_f} = [fold_sites{max_f}; uniq_site(sort_idx(s))];
        used_site(sort_idx(s)) = 1;
        sz_dif = N_fold - fold_sz;
    end
end

fold_subj = cell(Nfolds, 1);
for f = 1:Nfolds
    uniq_idx = [];
    for i = 1:length(fold_sites{f})
        uniq_idx = [uniq_idx; find(uniq_site == fold_sites{f}(i))];
    end
    
    fold_subj{f} = cat(1, site_subj{uniq_idx});
    fold_subj{f} = intersect(subjects, fold_subj{f}, 'stable');
end


end

