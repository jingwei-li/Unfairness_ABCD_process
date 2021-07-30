function ABCD_resample_split(outdir, subj_ls, orig_split_dir, csvname, bhvr_ls)

% ABCD_resample_split(outdir, subj_ls, orig_split_dir, csvname, bhvr_ls)
%
% Split the resampled subjects into folds.

subj_hdr = 'subjectkey';
site_hdr = 'site';
race_hdr = 'race';

if(~exist(outdir)); mkdir(outdir); end;

ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');
if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end

[subjects, nsub] = CBIG_text2cell(subj_ls);
[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
d = readtable(csvname);

all_idx = [];
for s = 1:nsub
    all_idx = [all_idx; find(strcmp(d.(subj_hdr), subjects{s}))];
end
AA = subjects(d.(race_hdr)(all_idx) == 2)';
WA = subjects(d.(race_hdr)(all_idx) == 1)';

AAidx = []; WAidx = [];
for s = 1:length(AA)    % AA WA should have the same length
    AAidx = [AAidx; find(strcmp(d.(subj_hdr), AA{s}))];
    WAidx = [WAidx; find(strcmp(d.(subj_hdr), WA{s}))];
end
siteAA = d.(site_hdr)(AAidx);
siteWA = d.(site_hdr)(WAidx);

for b = 1:nbhvr
    fprintf('%s\n', bhvr_nm{b})
    orig = load(fullfile(orig_split_dir, ['sub_fold_pass_rs_pass_pheno_' bhvr_nm{b} '.mat']));
    Nfolds = length(orig.sub_fold);

    for f = 1:Nfolds
        [~, ~, idx] = intersect(orig.sub_fold(f).subject_list, d.(subj_hdr), 'stable');
        curr_race = d.(race_hdr)(idx);
        nonAAWA = orig.sub_fold(f).subject_list(curr_race ~= 1 & curr_race ~= 2);
        curr_site = d.(site_hdr)(idx);
        uq_site = unique(curr_site);

        curr_AAidx = logical(zeros(length(AA), 1));
        curr_WAidx = logical(zeros(length(WA), 1));
        for st = 1:length(uq_site)
            curr_AAidx = curr_AAidx | strcmp(siteAA, uq_site{st});
            curr_WAidx = curr_WAidx | strcmp(siteWA, uq_site{st});
        end
        selAA = AA(curr_AAidx);
        selWA = WA(curr_WAidx);

        all_sel = [nonAAWA; selAA; selWA];
        [~, sel_idx] = ismember(all_sel, subjects);
        sub_fold(f).subject_list = subjects(sort(sel_idx))';
        sub_fold(f).fold_index = zeros(nsub, 1);
        for s = 1:nsub
            if(any(strcmp(sub_fold(f).subject_list, subjects{s})))
                sub_fold(f).fold_index(s) = 1;
            end
        end
    end
    save(fullfile(outdir, [bhvr_nm{b} '.mat']), 'sub_fold')
end
    
end