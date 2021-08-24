function [ subjects, subjdirs, goodRuns, goodSubjects ] = vgdl_getSubjectsDirsAndRuns()

% Get the list of subjects, subject directories, and number of runs for the fMRI GLM code
% copied from exploration_getSubjectsDirsAndRuns from https://github.com/tomov/Exploration-Data-Analysis
%


% the participant id as entered in psychopy
subjects = 1:32;

% should be identical to the list of subjects in the csv file
% and in the same order
% this is a basic assumption for getGoodSubjects() to work
% we are listing them here explicitly as a sanity check for the csv file
%
%assert(mean(strcmp(subjects, unique(data.participant)')) == 1);

% the names of the CORRESPONDING directories from CBS central
subjdirs = {'200311_VGDL_001', '200312_VGDL_002', '200313_VGDL_003', '200313_VGDL_004', ...
            '200314_VGDL_005', '200314_VGDL_006', '200315_VGDL_007', '200315_VGDL_008', ...
            '210121_VGDL_009', '210122_VGDL_010', '210126_VGDL_011', '210407_VGDL_012', ...
            '210411_VGDL_013', '212104_VGDL_014', '210422_VGDL_015', '210425_VGDL_016', ...
            '210425_VGDL_017', '210503_VGDL_018', '210505_VGDL_019', '20210521_VGDL_020', ...
            '210608_VGDL_021', '210614_VGDL_022', '20210616_VGDL_023', '210702_VGDL_024', ...
            '210703_VGDL_025', '210703_VGDL_026', '210703_VGDL_027', '210703_VGDL_028', ...
            '210705_VGDL_029', '210705_VGDL_030', '210705_VGDL_031', '210705_VGDL_032'};


% assumes runs are always in order: 1,2,3,4,...
%nRuns = {8,8}; % runs per subject

% which runs to include/exclude for each subject
goodRuns = {logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1 0]), logical([1 1 1 1 1 1]), logical([1 1 1 1 0 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 0 0 0]), logical([1 1 1 1 0 1]), ...
            logical([1 1 1 1 1 1]), logical([1 1 1 0 1 1]), logical([1 1 1 0 1 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([0 1 1 1 1 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1])};

% optionally, only use odd runs
% see GLM 11
%for s = 1:numel(goodRuns)
%    goodRuns{s}(2:2:end) = 0;
%end

% optionally, only use even runs
% see GLM 35
%for s = 1:numel(goodRuns)
%    goodRuns{s}(1:2:end) = 0;
%end


% which subjects are good
%
% bad subjects:
% 14 -- ambidextrous
% 17 -- brain abnormality in left thalamus
% 22 -- Prashant, lots of movement
%
% excluded runs:
% 11, run 5 -- game crashed after 111 TRs
% 15, runs 4,5,6 -- excessive motion
% 16, run 5 -- timeout popup for ~20 s
% 18, run 3 -- no anterior head coil after break
% 19, run 4 -- no anterior head coil after break
% 23, run 1 -- excessive motion
goodSubjects = [1:13 15:16 18:21 23:32];
 
assert(numel(subjects) == numel(subjdirs));
assert(numel(subjects) == numel(goodRuns));
