% learning timecourses based on interactions with approach/avoid/neutral objects

close all
clear all

EXPT=vgdl_expt;

conn = mongo('holy7c22109.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54')


frequency = 20; % Hz; IMPORTANT!!! must match frame_rate in vgdl/core.py
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
            learning_offset(g, subj_id).(valences{v}) = struct;
            learning_offset_collapsed(g, subj_id).(valences{v}) = [];
            learning_avn(g).(valences{v}) = [];
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
                for s=1:numel(sprite_names)
                    effects = avn.(valences{v}).(sprite_names{s})(avn.game_ids == g);
                    if any(isnan(effects))
                        % sprite from a different game => continue
                        assert(all(isnan(effects)), 'Either all or none of the effects should be NaN');
                        continue;
                    end
                    if ~ismember(sprite_names{s}, fieldnames(learning(g,subj_id).(valences{v})))
                        learning(g,subj_id).(valences{v}).(sprite_names{s}) = [];
                    end
                    learning(g,subj_id).(valences{v}).(sprite_names{s}) = [learning(g,subj_id).(valences{v}).(sprite_names{s}); effects];
                end

            end
        end

    end
end

% 2) offset timeseries from first event, for each game, subject, valance, sprite
for g = 1:6
    for subj_id = 1:length(EXPT.subject)
        for v=1:numel(valences)
            sprite_names=fieldnames(learning(g, subj_id).(valences{v}));

            % offset timecourses from first event
            max_length = 0;
            for s=1:numel(sprite_names)
                effects = learning(g, subj_id).(valences{v}).(sprite_names{s});
                offset = min(find(effects));
                if ~isempty(offset)
                    % offset effects, if there are any
                    effects = effects(offset:end);
                end
                learning_offset(g, subj_id).(valences{v}).(sprite_names{s}) = effects;
                max_length = max(max_length, length(learning_offset(g, subj_id).(valences{v}).(sprite_names{s})));
            end

            % pad all sprites to the same length
            for s=1:numel(sprite_names)
                remaining_length = max_length - length(learning_offset(g, subj_id).(valences{v}).(sprite_names{s}));
                learning_offset(g, subj_id).(valences{v}).(sprite_names{s}) = [learning_offset(g, subj_id).(valences{v}).(sprite_names{s}); nan(remaining_length, 1)];
            end
        end
    end
end

% 3) collapse across sprites and pad to the same length for all subjects
%    also smooth
for g = 1:6
    for v=1:numel(valences)
        % collate sprites and average
        max_length = 0;
        for subj_id = 1:length(EXPT.subject)
            sprite_names=fieldnames(learning(g, subj_id).(valences{v}));

            collated = [];
            for s=1:numel(sprite_names)
                effects = learning_offset(g, subj_id).(valences{v}).(sprite_names{s});
                collated = [collated, effects];
            end
            learning_offset_collapsed(g, subj_id).(valences{v}) = mean(collated, 2);
            %learning_offset_collapsed(g, subj_id).(valences{v}) = nanmean(collated, 2); % might have too much variance at the tails

            max_length = max(max_length, length(learning_offset_collapsed(g, subj_id).(valences{v})));
        end

        % pad for the same length for all subjects
        for subj_id = 1:length(EXPT.subject)
            remaining_length = max_length - length(learning_offset_collapsed(g, subj_id).(valences{v}));
            learning_offset_collapsed(g, subj_id).(valences{v}) = [learning_offset_collapsed(g, subj_id).(valences{v}); nan(remaining_length, 1)];
        end

        % smooth
        for subj_id = 1:length(EXPT.subject)
            learning_offset_collapsed_smooth(g, subj_id).(valences{v}) = imgaussfilt(learning_offset_collapsed(g, subj_id).(valences{v}), sigma);
        end
    end
end

% 4) collate across subjects
for g = 1:6
    for v=1:numel(valences)
        learning_avn(g).(valences{v}) = [learning_offset_collapsed(g,:).(valences{v})];
        learning_avn_smooth(g).(valences{v}) = [learning_offset_collapsed_smooth(g,:).(valences{v})];
    end
end


filename

clear conn;
save(filename);
