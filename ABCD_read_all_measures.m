function ABCD_read_all_measures(fmri_dir, subj_list, outdir, outstem)

dohist = 0;
subjectkey = CBIG_text2cell(subj_list);

%% read race
[race, race_hdr] = ABCD_read_race(subj_list, dohist);
race_empty = cellfun(@isempty, race);
fprintf('Empty race: %d subjects.\n', length(find(race_empty)))

%% read confounds
% 1. age
[age, age_hdr] = ABCD_read_age(subj_list, race, dohist);
age_empty = cellfun(@isempty, age);
fprintf('Empty age: %d subjects.\n', length(find(age_empty)));

% 2. sex
[sex, sex_hdr] = ABCD_read_sex(subj_list, race, dohist);
sex_empty = cellfun(@isempty, sex);
fprintf('Empty sex: %d subjects.\n', length(find(sex_empty)));

% 3. site
[site, site_hdr] = ABCD_read_site(subj_list, race, dohist);
site_empty = cellfun(@isempty, site);
fprintf('Empty site: %d subjects.\n', length(find(site_empty)));

% 4. motion
[FD, DVARS] = ABCD_read_motion(fmri_dir, subj_list, 1, outdir, outstem);

% 5. brain volume
[ICV] = ABCD_read_ICV(subj_list, race, dohist);
ICV_empty = isnan(ICV);
fprintf('Empty ICV: %d subjects.\n', length(find(ICV_empty)));

% 6. parental education
[peduc, peduc_comb, peduc_avg, peduc_hdr, peduc_colloquial] = ABCD_read_Prt_Educ(subj_list, race, dohist);
peduc_empty = isnan(peduc_avg);
fprintf('Empty parental education: %d subjects.\n', length(find(peduc_empty)));

subjectkey = subjectkey';
d = table(subjectkey, race, age, sex, site, FD, DVARS, ICV, peduc_avg);
for i = 1:length(peduc_hdr)
    d.(peduc_hdr{i}) = peduc(:,i);
end


%% read behaviors
behaviors = [];
colloquial = [];

%1. Rey Auditory Verval Learning Test
[RAVLT, RAVLT_hdr, RAVLT_colloquial] = ABCD_read_RAVLT(subj_list, race, dohist);
RAVLT_empty = isnan(sum(RAVLT,2));
fprintf('Empty RAVLT: %d subjects.\n', length(find(RAVLT_empty)));
for i = 1:length(RAVLT_hdr)
    d.(RAVLT_hdr{i}) = RAVLT(:,i);
end
behaviors = [behaviors RAVLT_hdr];
colloquial = [colloquial RAVLT_colloquial];

% 2. Wechsler Intelligence Scale for Children-V: Matrix reasoning
[WISC, WISC_hdr, WISC_colloquial] = ABCD_read_WISC(subj_list, race, dohist);
WISC_empty = isnan(sum(WISC,2));
fprintf('Empty WISC: %d subjects.\n', length(find(WISC_empty)));
for i = 1:length(WISC_hdr)
    d.(WISC_hdr{i}) = WISC(:,i);
end
behaviors = [behaviors WISC_hdr];
colloquial = [colloquial WISC_colloquial];

% 3. NIH Toolbox
[NIH, NIH_hdr, NIH_colloquial] = ABCD_read_NIH(subj_list, race, dohist);
NIH_empty = isnan(sum(NIH,2));
fprintf('Empty NIH: %d subjects.\n', length(find(NIH_empty)));
for i = 1:length(NIH_hdr)
    d.(NIH_hdr{i}) = NIH(:,i);
end
behaviors = [behaviors NIH_hdr];
colloquial = [colloquial NIH_colloquial];

% 4. Little Man Task
[LMT, LMT_hdr, LMT_colloquial] = ABCD_read_LittleMan(subj_list, race, dohist);
LMT_empty = isnan(sum(LMT,2));
fprintf('Empty LMT: %d subjects.\n', length(find(LMT_empty)));
for i = 1:length(LMT_hdr)
    d.(LMT_hdr{i}) = LMT(:,i);
end
behaviors = [behaviors LMT_hdr];
colloquial = [colloquial LMT_colloquial];

% 5. Aachenbach Child Behavior Check List
[CBCL, CBCL_hdr, CBCL_colloquial] = ABCD_read_CBCL(subj_list, race, dohist);
CBCL_empty = isnan(sum(CBCL,2));
fprintf('Empty CBCL: %d subjects.\n', length(find(CBCL_empty)));
for i = 1:length(CBCL_hdr)
    d.(CBCL_hdr{i}) = CBCL(:,i);
end
behaviors = [behaviors CBCL_hdr];
colloquial = [colloquial CBCL_colloquial];

% 6. Parental General Behavior Inventory: Mania
[PGBI, PGBI_hdr, PGBI_colloquial] = ABCD_read_PGBI(subj_list, race, dohist);
PGBI_empty = isnan(sum(PGBI,2));
fprintf('Empty PGBI: %d subjects.\n', length(find(PGBI_empty)));
for i = 1:length(PGBI_hdr)
    d.(PGBI_hdr{i}) = PGBI(:,i);
end
behaviors = [behaviors PGBI_hdr];
colloquial = [colloquial PGBI_colloquial];

% 7. Pediatric Psychosis Questionnaire
[PPS, PPS_hdr, PPS_colloquial] = ABCD_read_PPS(subj_list, race, dohist);
PPS_empty = isnan(sum(PPS,2));
fprintf('Empty PPS: %d subjects.\n', length(find(PPS_empty)));
for i = 1:length(PPS_hdr)
    d.(PPS_hdr{i}) = PPS(:,i);
end
behaviors = [behaviors PPS_hdr];
colloquial = [colloquial PPS_colloquial];

% 8. Modified UPPS-P for Children from PhenX
[UPPS, UPPS_hdr, UPPS_colloquial] = ABCD_read_UPPS(subj_list, race, dohist);
UPPS_empty = isnan(sum(UPPS,2));
fprintf('Empty UPPS: %d subjects.\n', length(find(UPPS_empty)));
for i = 1:length(UPPS_hdr)
    d.(UPPS_hdr{i}) = UPPS(:,i);
end
behaviors = [behaviors UPPS_hdr];
colloquial = [colloquial UPPS_hdr];

% 9. Behavioral Inhibition & Activation 
[BISBAS, BISBAS_hdr, BISBAS_colloquial] = ABCD_read_BISBAS(subj_list, race, dohist);
BISBAS_empty = isnan(sum(BISBAS,2));
fprintf('Empty BISBAS: %d subjects.\n', length(find(BISBAS_empty)));
for i = 1:length(BISBAS_hdr)
    d.(BISBAS_hdr{i}) = BISBAS(:,i);
end
behaviors = [behaviors BISBAS_hdr];
colloquial = [colloquial BISBAS_colloquial];


%% integrate results
any_empty = race_empty | age_empty | sex_empty | site_empty | ICV_empty | peduc_empty | ...
    RAVLT_empty | WISC_empty | NIH_empty | LMT_empty | CBCL_empty | PGBI_empty | ...
    PPS_empty | UPPS_empty | BISBAS_empty;
pass_subj = subjectkey(~any_empty);
CBIG_cell2text(pass_subj, fullfile(outdir, ['subjects' outstem '_pass_pheno.txt']));
writetable(d, fullfile(outdir, ['phenotypes' outstem '.txt']));
CBIG_cell2text(behaviors, fullfile(outdir, ['behavior_list.txt']))
CBIG_cell2text(colloquial, fullfile(outdir, ['colloquial_list.txt']))

end

