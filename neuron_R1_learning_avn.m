close all
clear all

EXPT=vgdl_expt;

conn = mongo('holy7c22109.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54')


frequency = 40; % Hz
level_duration = 60;  % s
sigma = 1 * frequency; % filter std

filename = fullfile(get_mat_dir(false), sprintf('neuron_R1_learning_avn_sigma=%.0f.mat', sigma));
filename

%lengths = [];
valences = {'approach','avoid','neutral'};

% Initialize
clear learning;
for g = 1:6
    for subj_id = 1:length(EXPT.subject)
        for v=1:numel(valences)
            learning(g, subj_id).(valences{v}) = struct;
            lengths(g, subj_id).(valences{v}) = 0;
        end
    end
end

% 1) get time series for each game, subject, valence, sprite
EXPT = vgdl_expt;
for subj_id = 1:length(EXPT.subject)
    for SPM_run_id = 1:length(EXPT.subject(subj_id).functional)

        run_id = get_behavioral_run_id(subj_id, SPM_run_id)
        fprintf('--- subj_id %d, SPM_run_id (behav %d) %d\n', subj_id, SPM_run_id, run_id);

        query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id) % in python we index runs from 0 (but not subjects) 
        run = find(conn, 'runs', 'query', query);

        avn = get_avn(subj_id, run, conn, true);

        game_ids = unique(avn.game_ids)';
        for g = game_ids

            for v=1:numel(valences)
                sprite_names=fieldnames(avn.(valences{v}));

                % appends sprites of corresponding valence for given run
                delta_length = nan;
                for s=1:numel(sprite_names)
                    if ~ismember(sprite_names{s}, fieldnames(learning(g,subj_id).(valences{v})))
                        learning(g,subj_id).(valences{v}).(sprite_names{s}) = nan(lengths(g,subj_id).(valences{v}), 1);
                    end
                    effects = avn.(valences{v}).(sprite_names{s});
                    delta_length = length(effects);
                    learning(g,subj_id).(valences{v}).(sprite_names{s}) = [learning(g,subj_id).(valences{v}).(sprite_names{s}); effects];
                end
                lengths(g, subj_id).(valences{v}) = lengths(g, subj_id).(valences{v}) + delta_length;

                % pad all sprites from other games with nan
                all_sprite_names = fieldnames(learning(g,subj_id).(valences{v}));
                remaining_sprite_names = all_sprite_names(~ismember(all_sprite_names, sprite_names));
                remaining_sprite_names
                delta_length
                for s=1:numel(remaining_sprite_names)
                    learning(g,subj_id).(valences{v}).(remaining_sprite_names{s}) = [learning(g,subj_id).(valences{v}).(remaining_sprite_names{s}); nan(delta_length,1)];
                end
            end
        end

        keyboard

    end
end

% 2) offset timeseries from first event

filename

clear conn;
save(filename);
