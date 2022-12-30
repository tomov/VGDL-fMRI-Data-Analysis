function [play_scores, play_wins, play_steps, play_durations, game_names, levels] = get_play_scores(conn, subj_id, run_ids, do_cache)
    % copy of get_instance_scores.m but for playsssss/episodes

    if ~exist('do_cache', 'var')
        do_cache = false;
    end
 
    filename = fullfile(get_mat_dir(false), sprintf('get_play_scores_subj%d_runs%s.mat', subj_id, sprintf('%d_', run_ids)));
    filename
    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    play_scores = [];
    play_wins = [];
    play_durations = [];
    play_steps = [];
    game_names = {};
    levels = [];

    for r = 1:length(run_ids)
        run_id = run_ids(r);

        q = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id);
        runs = find(conn, 'runs', 'query', q)
        if length(runs) == 0
            fprintf('no run %d for subject %d\n', run_id, subj_id);
            continue
        end
        assert(length(runs) == 1);
        run = runs(1);

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
                    proj = '{"win": 1.0, "score": 1.0, "end_time": 1.0, "start_time": 1.0, "actions": 1.0}';
                    plays = find(conn, 'plays', 'query', q, 'projection', proj);
                    assert(length(plays) == 1);
                    play = plays(1);

                    play_wins = [play_wins, play.win];
                    play_scores = [play_scores, play.score];
                    play_durations = [play_durations, play.end_time - play.start_time];
                    play_steps = [play_steps, length(play.actions)];
                    game_names = [game_names, {game_name}];
                    levels = [levels, instance.level_id + 1];
                end

            end
        end
    end

    if do_cache
        save(filename, 'play_scores', 'game_names', 'play_wins', 'play_durations', 'levels', '-v7.3');
    end
