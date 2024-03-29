function [subj_philips, isPhilips] = ABCD_get_subj_Philips(subj_list)

% [subj_philips, isPhilips] = ABCD_get_subj_Philips(subj_list)
%
% Subjects collected on Philips scanners were problematic:
% https://github.com/ABCD-STUDY/fMRI-cleanup.
% These subjects need to be excluded from the study.
% This function checks if any subject in the given list was collected on
% Philips scanners.
%
% Inputs:
% - subj_list
%   List of subjects with all required phenotypes, and passed rs-fMRI quality 
%   control (full path).
%
% Outputs:
% - subj_philips
%   The subject IDs which were collected on Philips scanners, if any.
%
% - isPhilips
%   A #subjects x 1 boolean vector. An entry of True indicates that this 
%   subject was collected on Philips scanners.
%
% Author: Jingwei Li 

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

mri_csv = '/mnt/isilon/CSC2/Yeolab/Data/ABCD/raw/documents/release2.0/ABCDstudyNDA/abcd_mri01.txt';
subj_hdr = 'subjectkey';
scanner_hdr = 'mri_info_manufacturer';
event_hdr = 'eventname';

d = readtable(mri_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');
scanner = cell(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        tmp_idx = tmp_idx & base_event;
        scanner(s) = d.(scanner_hdr)(tmp_idx);
    end
end

isPhilips = strcmp(scanner, 'Philips Medical Systems');
subj_philips = subjects(isPhilips);


end

