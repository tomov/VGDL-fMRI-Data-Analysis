function [legacy_fields, visuals, fields] = get_visuals(subj_id, run, conn, do_cache, plays_post_collection)

%clear all;
%conn = mongo('holy7c22105.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54', 'UserName', 'reader', 'Password', 'parolatamadafaqa')
%subj_id = 1;
%run_id = 6;
%query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id) % in python we index runs from 0 (but not subjects) 
%run = find(conn, 'runs', 'query', query)
%assert(length(run) == 1);
%do_cache = false;

    % helper function to get visual regressors for each frame in vgdl_create_multi
    % copied & modified from get_keypresses
    % note run is a struct
    %

    if ~exist('do_cache', 'var')
        do_cache = false;
    end

    % optionally cache
    if ~exist('plays_post_collection', 'var')
        plays_post_collection = 'plays_post';
        filename = fullfile(get_mat_dir(false), sprintf('get_visuals_new_subj%d_run%d.mat', subj_id, run.run_id));
    else
        filename = fullfile(get_mat_dir(false), sprintf('get_visuals_new_subj%d_run%d_c=%s.mat', subj_id, run.run_id, plays_post_collection));
    end
    filename

    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    % TODO tight coupling with vgdl_create_multi, case 10-20
    legacy_fields = {'timestamps', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', 'avatar_collision_flag', 'effectsByCol'};
    manual_fields = {'dscore', 'win', 'timeout', 'loss', 'durations', 'ended', 'play_start', 'play_end'};
    fields = [legacy_fields, manual_fields, {'score'}];

    visuals = struct;
    for i = 1:numel(fields)
        visuals.(fields{i}) = [];
    end

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

                plays_post = find(conn, plays_post_collection, 'query', q);
                assert(length(plays_post) == 1);
                play_post = plays_post(1);

                for i = 1:numel(fields)
                    if ~ismember(fields{i}, manual_fields)
                        visuals.(fields{i}) = [visuals.(fields{i}); play_post.(fields{i})(2:end-1,:)]; % skip first and last frame to align with output of get_regressors, b/c EMPA does too
                    end
                end

                % Manual fields
                %durations = play_post.timestamps(2:end) - play_post.timestamps(1:end-1);
                %durations = [durations; mean(durations)]; % guesstimate duration of last frame
                durations = play_post.timestamps(3:end) - play_post.timestamps(2:end-1); % skip first and last frame
                visuals.durations = [visuals.durations; durations];

                dscore = play_post.score(3:end) - play_post.score(2:end-1); % skip first and last frame
                visuals.dscore = [visuals.dscore; dscore];

                win = logical(zeros(size(play_post.win(2:end-1,:))));
                loss = logical(zeros(size(play_post.win(2:end-1,:))));
                timeout = logical(zeros(size(play_post.win(2:end-1,:))));
                ended = logical(zeros(size(play_post.win(2:end-1,:))));
                play_start = logical(zeros(size(play_post.win(2:end-1,:))));
                play_end = logical(zeros(size(play_post.win(2:end-1,:))));
                % Notice that here we actually take the last timestamp instead of the second to last one
                if length(win) > 0
                    if ~iscell(play_post.win) 
                        win_final = play_post.win(end);
                    else
                        win_final = play_post.win{end};
                    end
                    switch win_final
                        case 1
                            win(end) = 1;
                            loss(end) = 0;
                            timeout(end) = 0;
                        case 0
                            win(end) = 0;
                            loss(end) = 1;
                            timeout(end) = 0;
                        case -1
                            win(end) = 0;
                            loss(end) = 0;
                            timeout(end) = 1;
                        otherwise
                            % End of level
                            win(end) = 0;
                            loss(end) = 0;
                            timeout(end) = 1;
                            %assert(false, 'Invalid value for win_final');
                            disp('Invalid value for win_final');
                    end
                    ended(end) = 1;
                    play_start(1) = 1;
                    play_end(end) = 1;
                    visuals.win = [visuals.win; win];
                    visuals.loss = [visuals.loss; loss];
                    visuals.timeout = [visuals.timeout; timeout];
                    visuals.ended = [visuals.ended; ended];
                    visuals.play_start = [visuals.play_start; play_start];
                    visuals.play_end = [visuals.play_end; play_end];
                end
            end
        end

    end

    visuals.timestamps = visuals.timestamps - run.scan_start_ts;


    if do_cache
        save(filename, 'visuals', 'legacy_fields', 'fields', '-v7.3');
    end
