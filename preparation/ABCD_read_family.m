function [fam_id, fam_mem] = ABCD_read_family(subj_list)

addpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))


fam_csv = '/mnt/eql/yeo12/data/ABCD/documents/release2.0/ABCDstudyNDA/acspsw03.txt';
fam_hdr = 'rel_family_id';
subj_hdr = 'subjectkey';
event_hdr = 'eventname';

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt';
end
[subjects, nsub] = CBIG_text2cell(subj_list);

% format in subj_list: e.g. NDARINV007W6H7B
% format in csv: e.g. "NDAR_INV007W6H7B"
subjects_csv = cell(nsub, 1);
for s = 1:nsub
    subjects_csv{s} = [subjects{s}(1:4) '_' subjects{s}(5:end)];
end

d = readtable(fam_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');
fam_id_read = d.(fam_hdr);

% select only the rows corresponding to required subjects
fam_id = cell(nsub, 1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        tmp_idx = tmp_idx & base_event;
        fam_id(s) = fam_id_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, fam_id);
fam_id(empty_idx) = {'NaN'};
fam_id = cellfun(@str2num, fam_id);

uniq_fam = unique(fam_id);
fam_mem = cell(length(uniq_fam), 1);
fam_size = zeros(length(uniq_fam), 1);
for f = 1:length(uniq_fam)
    fam_mem{f} = subjects(fam_id == uniq_fam(f));
    fam_size(f) = length(fam_mem{f});
end

rmpath(genpath( '/data/users/jingweil/storage/from_HOME/code/plotting_functions/'))

end

