function [instance_scores, instance_wins, instance_success, instance_success_rates, game_names, levels] = get_instance_scores(conn, subj_id, run_ids, do_cache)

    % get score for each instance (i.e. level) as max across all won plays 
    % this was the same way we determined to pay out for the fMRI study

    if ~exist('do_cache', 'var')
        do_cache = false;
    end
 
    filename = fullfile(get_mat_dir(false), sprintf('get_instance_scores_subj%d_runs%s.mat', subj_id, sprintf('%d_', run_ids)));
    filename
    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    instance_scores = [];
    instance_wins = [];
    instance_success = [];
    instance_success_rates = [];
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

                instance_score = 0;
                instance_win = 0;

                for p = 1:nplays
                    % fetch plays one by one, b/c o/w OOM
                    q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d, "play_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id, p - 1);
                    proj = '{"win": 1.0, "score": 1.0}';
                    plays = find(conn, 'plays', 'query', q, 'projection', proj);
                    assert(length(plays) == 1);
                    play = plays(1);

                    if play.win
                        instance_score = max(instance_score, play.score);
                        instance_win = instance_win + 1;
                    end
                end

                instance_scores = [instance_scores, instance_score];
                instance_wins = [instance_wins, instance_win];
                instance_success = [instance_success, instance_win > 0];
                instance_success_rates = [instance_success_rates, instance_win / nplays];
                game_names = [game_names, {game_name}];
                levels = [levels, instance.level_id + 1];
            end
        end
    end

    if do_cache
        save(filename, 'instance_scores', 'game_names', 'instance_wins', 'instance_success', 'instance_success_rates', 'levels', '-v7.3');
    end
