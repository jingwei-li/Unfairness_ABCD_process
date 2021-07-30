function ABCD_upsampleAA_downsampleWA(outdir, subj_ls, csvname )

% ABCD_upsampleAA_downsampleWA(outdir, subj_ls, csvname )
%
% 

race_hdr = 'race';
subj_hdr = 'subjectkey';
niter = 100;

ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');
if(~exist('subj_ls', 'var') || isempty(subj_ls))
    subj_ls = fullfile(ls_dir, 'subjects_pass_rs_pass_pheno.txt');
end
if(~exist('csvname', 'var') || isempty(csvname))
    csvname = fullfile(ls_dir, 'phenotypes_pass_rs.txt');
end

[subjects, nsub] = CBIG_text2cell(subj_ls);
d = readtable(csvname);
[~, ~, idx] = intersect(subjects, d.(subj_hdr), 'stable');

race = d.(race_hdr)(idx);
AA = subjects(race==2);
WA = subjects(race==1);
target = (length(AA) + length(WA)) / 2;
AAWAidx = find(race==2 | race == 1);

if(~exist(outdir)); mkdir(outdir); end;
rng('default'); rng(1000);
for i = 1:niter
    spAAidx = sort(datasample(1:length(AA), target, 'replace', true), 'ascend');
    spWAidx = sort(datasample(1:length(WA), target, 'replace', false), 'ascend');

    spAA = AA(spAAidx);
    spWA = WA(spWAidx);
    new_subj = subjects;
    new_subj(AAWAidx(1:target)) = spAA;
    new_subj(AAWAidx(target+1:end)) = spWA;

    CBIG_cell2text(new_subj, fullfile(outdir, ['subjects_sample' num2str(i) '.txt']))
end
    
end