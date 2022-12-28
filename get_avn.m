function [avn, fields] = get_avn(subj_id, run, conn, do_cache, approach_avoid_collection)

%clear all;
%conn = mongo('holy7c22109.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54')
%subj_id = 1;
%run_id = 1;
%query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id) % in python we index runs from 0 (but not subjects) 
%run = find(conn, 'runs', 'query', query)
%assert(length(run) == 1);
%do_cache = false;

    % get approach/avoid; goes together with fmri_approachAvoid.py

    if ~exist('do_cache', 'var')
        do_cache = false;
    end

    % optionally cache
    if ~exist('approach_avoid_collection', 'var')
        approach_avoid_collection = 'approach_avoid';
        filename = fullfile(get_mat_dir(false), sprintf('get_avn_subj%d_run%d.mat', subj_id, run.run_id));
    else
        filename = fullfile(get_mat_dir(false), sprintf('get_avn_subj%d_run%d_c=%s.mat', subj_id, run.run_id, approach_avoid_collection));
    end
    filename

    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end


    valences = {'approach','avoid','neutral'};
    fields = valences;


    game_names_ordered = get_game_names_ordered(subj_id);
    game_name_to_id = containers.Map(game_names_ordered, 1:6);

    avn = struct;
    for i = 1:numel(valences)
        avn.(valences{i}) = struct;
    end
    avn.game_ids = [];
    avn.block_ids = [];
    avn.instance_ids = [];
    avn.play_ids = [];
    avn.state_timestamps = []; % these all come directly from states, so they are state['ts'] timestamps;

    avn.game_names_ordered = game_names_ordered;
    avn.game_name_to_id = game_name_to_id;

    blocks = run.blocks;
    for b = 1:length(blocks)
        block = blocks(b);
        instances = block.instances;

        game_name = block.game.name;
        
        for i = 1:length(instances)
            instance = instances(i);

            q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id);
            nplays = count(conn, 'plays', 'query', q);

            for p = 1:nplays
                % fetch plays one by one, b/c o/w OOM
                q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d, "play_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id, p - 1);
                plays = find(conn, 'plays', 'query', q);
                assert(length(plays) == 1);
                play = plays(1);

                approach_avoids = find(conn, approach_avoid_collection, 'query', q);
                assert(length(approach_avoids) == 1);
                approach_avoid = approach_avoids(1);

                state_timestamps = approach_avoid.state_timestamps(2:end-1);

                for v=1:numel(valences)
                    sprite_names = fieldnames(approach_avoid.effects_by_valence.(valences{v}));

                    % appends sprites of corresponding valence for given play
                    for s=1:numel(sprite_names)
                        if ~ismember(sprite_names{s}, fieldnames(avn.(valences{v})))
                            % new sprite type => pad with nans
                            avn.(valences{v}).(sprite_names{s}) = nan(size(avn.state_timestamps));
                        end
                        effects = approach_avoid.effects_by_valence.(valences{v}).(sprite_names{s})(2:end-1);
                        avn.(valences{v}).(sprite_names{s}) = [avn.(valences{v}).(sprite_names{s}); effects];
                    end

                    % pad all sprites from other games with nan
                    all_sprite_names = fieldnames(avn.(valences{v}));
                    remaining_sprite_names = all_sprite_names(~ismember(all_sprite_names, sprite_names));
                    for s=1:numel(remaining_sprite_names)
                        avn.(valences{v}).(remaining_sprite_names{s}) = [avn.(valences{v}).(remaining_sprite_names{s}); nan(size(state_timestamps))];
                    end
                end

                avn.state_timestamps = [avn.state_timestamps; state_timestamps]; % it's important that this happens after the sprite updates
                avn.game_ids = [avn.game_ids; ones(size(state_timestamps)) * game_name_to_id(game_name)]; 
                avn.block_ids = [avn.block_ids; ones(size(state_timestamps)) * block.block_id]; 
                avn.instance_ids = [avn.instance_ids; ones(size(state_timestamps)) * instance.instance_id];
                avn.play_ids = [avn.play_ids; ones(size(state_timestamps)) * play.play_id];

            end
        end

    end

    avn.state_timestamps = avn.state_timestamps - run.scan_start_ts;



    if do_cache
        save(filename, 'avn', 'fields', '-v7.3');
    end
