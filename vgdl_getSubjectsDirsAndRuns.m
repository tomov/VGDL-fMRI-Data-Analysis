function [ subjects, subjdirs, goodRuns, goodSubjects ] = vgdl_getSubjectsDirsAndRuns()

% Get the list of subjects, subject directories, and number of runs for the fMRI GLM code
% copied from exploration_getSubjectsDirsAndRuns from https://github.com/tomov/Exploration-Data-Analysis
%


% the participant id as entered in psychopy
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];

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
            '210411_VGDL_013'};


% assumes runs are always in order: 1,2,3,4,...
%nRuns = {8,8}; % runs per subject

% which runs to include/exclude for each subject
goodRuns = {logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1]),   logical([1 1 1 1 1 1]), logical([1 1 1 1 0 1]), logical([1 1 1 1 1 1]), ...
            logical([1 1 1 1 1 1])};

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
goodSubjects = 1:13;
 
assert(numel(subjects) == numel(subjdirs));
assert(numel(subjects) == numel(goodRuns));
