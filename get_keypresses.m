function [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn)

    % helper function to get keypresses, keyholds, etc in vgdl_create_multi
    % note run is a struct
    %

    keyNames = [];

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

                if isempty(keyNames)
                    keyNames = fieldnames(play.keyholds);
                    keyholds = cell(1,length(keyNames));
                    keyholds_post = cell(1,length(keyNames));
                    keypresses = cell(1,length(keyNames));
                end

                % key hold boxcar regressors
                for k = 1:numel(keyNames)
                    kh = play.keyholds.(keyNames{k});
                    if length(kh) > 0
                        kh(:,1) = kh(:,1) - run.scan_start_ts;
                    end
                    keyholds{k} = [keyholds{k}; kh];

                    kh_post = play_post.keyholds.(keyNames{k});
                    if length(kh_post) > 0
                        kh_post(:,1) = kh_post(:,1) - run.scan_start_ts;
                    end
                    keyholds_post{k} = [keyholds_post{k}; kh_post];

                    kp = play_post.keypresses.(keyNames{k});
                    if length(kp) > 0
                        kp(:,1) = kp(:,1) - run.scan_start_ts;
                    end
                    keypresses{k} = [keypresses{k}; kp];
                end
            end
        end

    end
