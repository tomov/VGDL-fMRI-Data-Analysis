function [onoff] = get_onoff(subj_id, run, conn, do_cache)

    % helper function to get screen on/off regressors in vgdl_create_multi
    % copied & modified from get_visuals
    % note run is a struct
    %

    if ~exist('do_cache', 'var')
        do_cache = false;
    end

    % optionally cache
    filename = sprintf('mat/get_onoff_subj%d_run%d.mat', subj_id, run.run_id);
    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    onoff = struct;
    onoff.block_start = [];
    onoff.block_end = [];
    onoff.instance_start = [];
    onoff.instance_end = [];
    onoff.play_start = [];
    onoff.play_end = [];

    blocks = run.blocks;
    for b = 1:length(blocks)
        block = blocks(b);
        instances = block.instances;

        game_name = block.game.name;

        onoff.block_start = [onoff.block_start; block.start_time - run.scan_start_ts];
        onoff.block_end = [onoff.block_end; block.end_time - run.scan_start_ts];
        
        for i = 1:length(instances)
            instance = instances(i);

            onoff.instance_start = [onoff.instance_start; instance.start_time - run.scan_start_ts];
            onoff.instance_end = [onoff.instance_end; instance.end_time - run.scan_start_ts];
        
            q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id);
            nplays = count(conn, 'plays', 'query', q);

            for p = 1:nplays
                % fetch plays one by one, b/c o/w OOM
                q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d, "play_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id, p - 1);
                plays = find(conn, 'plays', 'query', q);
                assert(length(plays) == 1);
                play = plays(1);

                assert(play.run_start_ts == run.scan_start_ts);
                onoff.play_start = [onoff.play_start; play.start_time - run.scan_start_ts];
                onoff.play_end = [onoff.play_end; play.end_time - run.scan_start_ts];
        
            end
        end

    end


    if do_cache
        save(filename, 'onoff', '-v7.3');
    end
