function [theory_change_flag_onsets] = get_regressors(subj_id, run, conn)

    % helper function to get regressors for each frame in vgdl_create_multi
    % copied & modified from GLM 3
    % note run is a struct
    %

    theory_change_flag_onsets = [];

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
                % one by one, b/c o/w OOM
                q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d, "play_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id, p - 1);
                plays = find(conn, 'plays', 'query', q);
                assert(length(plays) == 1);
                play = plays(1);

                regressors = find(conn, 'regressors', 'query', q);
                %regressors = find(conn, 'regressors', 'query', q, 'sort', '{"dt": -1.0}'); % momchil: assume latest one is the correct one 
                if length(regressors) == 0
                    % TODO remove once we fix plaqueAttack & lemmings!
                    fprintf('       NO REGRESSORS FOR %s: %d %d %d %d\n', play.game_name, play.run_id, play.block_id, play.instance_id, play.play_id);
                    continue;
                end
                assert(length(regressors) == 1);
                reg = regressors(1);

                tc = reg.regressors.theory_change_flag;
                for i = 1:length(tc)
                    if tc{i}{1} % theory changed
                        theory_change_flag_onsets = [theory_change_flag_onsets tc{i}{3} - run.scan_start_ts];
                    end
                end
            end

        end

    end
