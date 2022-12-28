close all
clear all

EXPT=vgdl_expt;

conn = mongo('holy7c22109.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54')


frequency = 40; % Hz
level_duration = 60;  % s
sigma = 1 * frequency; % filter std

filename = fullfile(get_mat_dir(false), sprintf('neuron_R1_learning_sigma=%.0f.mat', sigma));
filename

%lengths = [];

% Initialize
clear learning;
for g = 1:6
    for l = 1:9
        learning(g, l).tcf_stretched = nan(length(EXPT.subject), frequency * level_duration);
        learning(g, l).tcf_smooth = nan(length(EXPT.subject), frequency * level_duration);
    end
end

%
EXPT = vgdl_expt;
for subj_id = 1:length(EXPT.subject)
    for SPM_run_id = 1:length(EXPT.subject(subj_id).functional)

        run_id = get_behavioral_run_id(subj_id, SPM_run_id)
        fprintf('--- subj_id %d, SPM_run_id (behav %d) %d\n', subj_id, SPM_run_id, run_id);

        query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id) % in python we index runs from 0 (but not subjects) 
        run = find(conn, 'runs', 'query', query);

        regs = get_regressors(subj_id, run, conn, true);

        game_ids = unique(regs.game_ids)';
        for g = game_ids
            instance_ids = unique(regs.instance_ids(regs.game_ids == g))';
            for i = instance_ids
                level = get_level_from_run_and_instance(run_id, i);
                fprintf('  ------------ game %d, level %d\n', g, level);

                tcf = regs.theory_change_flag((regs.game_ids == g) & (regs.instance_ids == i));

                %lengths = [lengths, length(tcf)];
                tcf_stretched = imresize(tcf, [frequency * level_duration, 1], 'nearest');
                tcf_smooth = imgaussfilt(tcf_stretched, sigma);
                
                learning(g, level).tcf_stretched(subj_id, :) = tcf_stretched;
                learning(g, level).tcf_smooth(subj_id, :) = tcf_smooth;
            end
        end
        %snaohtu
    end
end

filename

clear conn;
save(filename);
