close all
clear all

% neuron_R1_learning.m, except without splitting across games

[subjects, subjdirs, goodRuns, goodSubjects ] = vgdl_getSubjectsDirsAndRuns();

conn = mongo('holy7c22109.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54')


frequency = 40;  % Hz
level_duration = 60;  % s
num_runs = 6;
levels_per_run = 9;
run_duration = level_duration * levels_per_run;

sigma = 10 * frequency; % filter std
%lengths = [];

filename = fullfile(get_mat_dir(false), sprintf('neuron_R1_learning_2_sigma=%.0f.mat', sigma));
filename


EXPT = vgdl_expt;
for subj_id = 1:length(EXPT.subject)
    learning(subj_id).tcf_stretched = [];
    learning(subj_id).tcf_smooth = [];

    for run_id = 1:num_runs

        query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id) % in python we index runs from 0 (but not subjects) 
        run = find(conn, 'runs', 'query', query);

        if goodRuns{subj_id}(run_id)
            regs = get_regressors(subj_id, run, conn, true);

            tcf = regs.theory_change_flag;

            tcf_stretched = imresize(tcf, [frequency * run_duration, 1], 'nearest');
            tcf_smooth = imgaussfilt(tcf_stretched, sigma);
            
            learning(subj_id).tcf_stretched = [learning(subj_id).tcf_stretched; tcf_stretched];
            learning(subj_id).tcf_smooth = [learning(subj_id).tcf_smooth; tcf_smooth];
        else
            % missing/bad run
            learning(subj_id).tcf_stretched = [learning(subj_id).tcf_stretched; nan(frequency*run_duration, 1)];
            learning(subj_id).tcf_smooth = [learning(subj_id).tcf_smooth; nan(frequency*run_duration, 1)];
        end
    end
end

filename

clear conn;
save(filename);
