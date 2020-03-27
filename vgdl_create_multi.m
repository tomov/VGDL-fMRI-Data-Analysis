function multi = vgdl_create_multi(glmodel, subj_id, run_id, save_output)
%glmodel = 2;
%subj_id = 1;
%run_id = 1;

clear multi;
save_output = true;

    % Create multi structure, helper function for creating EXPT in
    % vgdl_expt.m
    % copied from exploration_create_multi.m from https://github.com/tomov/Exploration-Data-Analysis
    %
    % USAGE: multi = vgdl_create_multi(model,subj_id,run_id)
    %
    % INPUTS:
    %   glmodel - positive integer indicating general linear model
    %   subj_id - integer specifying which subject is being analyzed
    %   run_id - integer specifying the run_id (1-indexed)
    %

    % OUTPUTS:
    %   multi - a sctructure with the folloowing fields
    %        .names{i}
    %        .onsets{i}
    %        .duratlsions{i}
    %        optional:
    %        .pmod(i).name
    %        .pmod(i).param
    %        .pmod(i).poly
    %
    % Momchil Tomov, Mar 2020

    if ~exist('save_output', 'var')
        save_output = false;
    end

    fprintf('glm %d, subj_id %d, run_id %d\n', glmodel, subj_id, run_id);


    [allSubjects, subj_dirs, goodRuns, goodSubjs] = vgdl_getSubjectsDirsAndRuns();

    SPM_run_id = run_id; % save the SPM run_id (SPM doesn't see bad run_ids)
  
    % skip bad run_ids
    % TODO enable
    %run_ids = find(goodRuns{find(subj_id == allSubjects)});
    %run_id = run_ids(run_id);
    fprintf('run_id %d \n', run_id);

    filename = sprintf('vgdl_create_multi_glm%d_subj%d_run%d.mat', glmodel, subj_id, run_id)

    % hack to make it work on the cluster until they install MongoDB toolbox
    % just pre-generate the multi's locally, then load them on the cluster
    try
        conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')
    catch e
        e
        fprintf('loading from %s\n', filename);
        load(filename);
        return
    end


    query = sprintf('{"subj_id": "%d"}', subj_id)

    subj = find(conn, 'subjects', 'query', query)
    assert(length(subj) > 0, 'Subject not found -- wrong subj_id?');
    assert(length(subj) <= 1, 'Too many subjects -- duplicate entries? Yikes!');

    query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id) % in python we index runs from 0 (but not subjects) 

    run = find(conn, 'runs', 'query', query)
    assert(length(run) == 1);
    
    % TODO nuisance regressors
    % start/stop play/instance/block
    % key down / key ups
    % visual/state prediction error (# of changed pixels ~= # of moved sprites)
    % # of collisions

    % null RSA:
    % - pixels
    % - hidden units of VAE trained to predict next scene
    % - AlexNet layers
    % 

    % GLMs
    %
    switch glmodel

        % condition = game, boxcars over blocks
        % look for systematic differences across games for things like agency, etc
        % e.g. contrast: 'chase + helper - butterflies - aliens'
        %
        case 1 

            idx = 0;

            keys = [];

            blocks = run.blocks;
            for b = 1:length(blocks)
                block = blocks(b);
                instances = block.instances;

                game_name = block.game.name;
                
                onsets = [];
                durs = [];
                for i = 1:length(instances)
                    instance = instances(i);

                    st = instance.start_time - run.scan_start_ts;
                    en = instance.end_time - run.scan_start_ts; % includes "YOU WON / LOST / etc" screen

                    onsets = [onsets, st];
                    durs = [durs, en - st];
                end

                % instance boxcar regressor
                idx = idx + 1;
                multi.names{idx} = game_name;
                multi.onsets{idx} = onsets;
                multi.durations{idx} = durs;
            end


        % key holds vs. button presses vs. key holds from button presses   NOT A REAL GLM for viz only
        % nuisance regressors
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        % also to compare keyholds_post computed from keypresses with actual keyholds
        %
        case 2

            idx = 0;

            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn);

            for k = 1:numel(keyNames)
                if size(keyholds{k}, 1) > 0
                    idx = idx + 1;
                    multi.names{idx} = ['keyholds_', keyNames{k}];
                    multi.onsets{idx} = keyholds{k}(:,1)';
                    multi.durations{idx} = keyholds{k}(:,2)';
                end

                if size(keyholds_post{k}, 1) > 0
                    idx = idx + 1;
                    multi.names{idx} = ['keyholds_post_', keyNames{k}];
                    multi.onsets{idx} = keyholds_post{k}(:,1)';
                    multi.durations{idx} = keyholds_post{k}(:,2)';
                end

                if size(keypresses{k}, 1) > 0
                    idx = idx + 1;
                    multi.names{idx} = ['keypresses_', keyNames{k}];
                    multi.onsets{idx} = keypresses{k}(:,1)';
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));
                end
            end



        % theory_change_flag 
        %
        case 3

            idx = 0;

            onsets = [];
            durs = [];

            blocks = run.blocks;
            for b = 1:length(blocks)
                block = blocks(b);
                instances = block.instances;

                game_name = block.game.name;
                
                for i = 1:length(instances)
                    instance = instances(i);

                    q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id);
                    plays = find(conn, 'plays', 'query', q);

                    for p = 1:length(plays)
                        play = plays(p);

                        q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d, "play_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id, play.play_id);
                        % momchil: assume latest one is the correct one TODO cleanup
                        %regressors = find(conn, 'regressors', 'query', q, 'sort', '{"dt": -1.0}');

                        tc = regressors(1).regressors.theory_change_flag;
                        for i = 1:length(tc)
                            if tc{i}{1} % theory changed
                                onsets = [onsets tc{i}{3} - run.scan_start_ts];
                                durs = [durs 0];
                            end
                        end
                    end
                end

            end

            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = durs;



        % same as 1 but with nuisance regressors
        % condition = game, boxcars over blocks
        % look for systematic differences across games for things like agency, etc
        % e.g. contrast: 'chase + helper - butterflies - aliens'
        %
        case 4 

            idx = 0;

            keys = [];

            blocks = run.blocks;
            for b = 1:length(blocks)
                block = blocks(b);
                instances = block.instances;

                game_name = block.game.name;
                
                onsets = [];
                durs = [];
                for i = 1:length(instances)
                    instance = instances(i);

                    st = instance.start_time - run.scan_start_ts;
                    en = instance.end_time - run.scan_start_ts; % includes "YOU WON / LOST / etc" screen

                    q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id);
                    nplays = count(conn, 'plays', 'query', q);

                    for p = 1:nplays
                        % one by one, b/c o/w OOM
                        q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d, "play_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id, p - 1);
                        plays = find(conn, 'plays', 'query', q);
                        assert(length(plays) == 1);
                        play = plays(1);

                        % TODO continue
                    end

                    onsets = [onsets, st];
                    durs = [durs, en - st];
                end

                % instance boxcar regressor
                idx = idx + 1;
                multi.names{idx} = game_name;
                multi.onsets{idx} = onsets;
                multi.durations{idx} = durs;
            end


        % button hold boxcars, (copied) subset of GLM 2
        % nuisance regressors
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        %
        case 5

            idx = 0;

            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn);

            if subj_id == 1
                % we screwed up keyholds for subject 1, so we use estimates from keypresses
                keyholds = keyholds_post
            end

            % key hold boxcar regressors
            for k = 1:numel(keyNames)
                if size(keyholds{k}, 1) > 0
                    idx = idx + 1;
                    multi.names{idx} = keyNames{k};
                    multi.onsets{idx} = keyholds{k}(:,1)';
                    multi.durations{idx} = keyholds{k}(:,2)';
                end
            end


        % key presses, (copied) subset of GLM 2, compare to GLM 5
        % nuisance regressors
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        %
        case 6

            idx = 0;

            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn);

            for k = 1:numel(keyNames)
                if size(keypresses{k}, 1) > 0
                    idx = idx + 1;
                    multi.names{idx} = ['keypresses_', keyNames{k}];
                    multi.onsets{idx} = keypresses{k}(:,1)';
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));
                end
            end




        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

    if save_output
        save(filename, 'multi', '-v7.3'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
    end


