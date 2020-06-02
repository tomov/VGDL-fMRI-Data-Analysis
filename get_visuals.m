function [fields, visuals] = get_visuals(subj_id, run, conn, do_cache)

    % helper function to get visual regressors for each frame in vgdl_create_multi
    % copied & modified from get_keypresses
    % note run is a struct
    %

    if ~exist('do_cache', 'var')
        do_cache = false;
    end

    % optionally cache
    filename = sprintf('mat/get_visuals_subj%d_run%d.mat', subj_id, run.run_id);
    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    % TODO tight coupling with vgdl_create_multi, case 10-20
    fields = {'timestamps', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', 'avatar_collision_flag', 'effectsByCol'};
    visuals = struct;
    for i = 1:numel(fields)
        visuals.(fields{i}) = [];
    end
    visuals.durations = [];

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

                plays_post = find(conn, 'plays_post', 'query', q);
                assert(length(plays_post) == 1);
                play_post = plays_post(1);

                for i = 1:numel(fields)
                    visuals.(fields{i}) = [visuals.(fields{i}); play_post.(fields{i})];
                end

                durations = play_post.timestamps(2:end) - play_post.timestamps(1:end-1);
                durations = [durations; mean(durations)]; % guesstimate duration of last frame
                visuals.durations = [visuals.durations; durations];
            end
        end

    end

    visuals.timestamps = visuals.timestamps - run.scan_start_ts;


    if do_cache
        save(filename, 'visuals', 'fields', '-v7.3');
    end
