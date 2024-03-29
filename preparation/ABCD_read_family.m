function [fam_id, fam_mem, grp_id] = ABCD_read_family(subj_list)

% [fam_id, fam_mem, grp_id] = ABCD_read_family(subj_list)
%
% Read family IDs of each subject from ABCD csv file.
%
% Inputs:
% - subj_list
%   List of subjects which passed fMRI prepreocessing quality control (full path). Default:
%   '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt'
%
% Outputs:
% - fam_id
%   A #subjects x 1 vector. Each array is the family ID of a subject.
% 
% - fam_mem
%   A cell array whose length is the number of unique families. Each array contains the 
%   subject IDs belong to the current family.
%
% - grp_id
%   A #subjects x 1 cell. Each subject was given a group id by ABCD to specify the 
%   twins/sibling relationships with the other members of the same family. This variable
%   will be used for multi-level block permutation test later.
%
% Author: Jingwei Li

fam_csv = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/raw/documents/release2.0/ABCDstudyNDA/acspsw03.txt';
fam_hdr = 'rel_family_id';
grp_hdr = 'rel_group_id';
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
grp_id_read = d.(grp_hdr);

% select only the rows corresponding to required subjects
fam_id = cell(nsub, 1);
grp_id = cell(nsub, 1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        tmp_idx = tmp_idx & base_event;
        fam_id(s) = fam_id_read(tmp_idx,:);
        grp_id(s) = grp_id_read(tmp_idx,:);
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

end

