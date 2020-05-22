function multi = vgdl_create_multi(glmodel, subj_id, run_id, save_output)
%glmodel = 7;
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

    filename = sprintf('mat/vgdl_create_multi_glm%d_subj%d_run%d.mat', glmodel, subj_id, run_id)

    % hack to make it work on the cluster until they install MongoDB toolbox
    % just pre-generate the multi's locally, then load them on the cluster
    if ~ismember(glmodel, 23)
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
    end
    
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

            [game_names, onsets, durs] = get_games(subj_id, run, conn);

            for i = 1:numel(game_names)
                % instance boxcar regressor
                idx = idx + 1;
                multi.names{idx} = game_names{i};
                multi.onsets{idx} = onsets{i};
                multi.durations{idx} = durs{i};
            end


        % key holds vs. button presses vs. key holds from button presses   NOT A REAL GLM for viz only
        % nuisance regressors
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        % also to compare keyholds_post computed from keypresses with actual keyholds
        %
        case 2

            idx = 0;

            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn, true);

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

            regs = get_regressors(subj_id, run, conn);
            onsets = regs.theory_change_flag_onsets;

            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;



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

            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn, true);

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

            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn, true);

            for k = 1:numel(keyNames)
                if size(keypresses{k}, 1) > 0
                    idx = idx + 1;
                    multi.names{idx} = ['keypresses_', keyNames{k}]; % TODO rm prefix... but already ran
                    multi.onsets{idx} = keypresses{k}(:,1)';
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));
                end
            end


        % frame nuisance regressors: visual control regressors for each frame
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        %
        case 7

            [fields, visuals] = get_visuals(subj_id, run, conn, true);

            multi.names{1} = 'frames';
            multi.onsets{1} = visuals.timestamps';
            multi.durations{1} = visuals.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            idx = 0;

            for i = 1:numel(fields)
                if ~ismember(fields{i}, {'timestamps', 'durations'}) 
                    if all(visuals.(fields{i}) == visuals.(fields{i})(1))
                        % constant
                        continue
                    end
                    idx = idx + 1;
                    multi.pmod(1).name{idx} = fields{i};
                    multi.pmod(1).param{idx} = visuals.(fields{i});
                    multi.pmod(1).poly{idx} = 1;
                end
            end

        % on/off nuisance regressors: visual control regressors for start/end block, instance, play
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        %
        case 8

            [onoff] = get_onoff(subj_id, run, conn, true);
            fields = fieldnames(onoff);

            idx = 0;

            for i = 1:numel(fields)
                idx = idx + 1;
                multi.names{idx} = fields{i};
                multi.onsets{idx} = onoff.(fields{i});
                multi.durations{idx} = zeros(size(multi.onsets{idx}));
            end


        % condition = game, boxcars over blocks (GLM 1)
        % + nuisance regressors (buttons, frames, on/off) (GLMs 5, 7, 8)
        % look for systematic differences across games for things like agency, etc
        % e.g. contrast: 'chase + helper - butterflies - aliens'
        % 
        case 9

            idx = 0;

            % from GLM 1: game instance boxcar regressor
            %
            [game_names, onsets, durs] = get_games(subj_id, run, conn);
            for i = 1:numel(game_names)
                idx = idx + 1;
                multi.names{idx} = game_names{i};
                multi.onsets{idx} = onsets{i};
                multi.durations{idx} = durs{i};
            end

            % from GLM 5: keyholds nuisance regressors
            %
            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn, true);
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

            % GLM 7: frame nuisance regressors
            %
            [fields, visuals] = get_visuals(subj_id, run, conn, true);
            idx = idx + 1;
            multi.names{idx} = 'frames';
            multi.onsets{idx} = visuals.timestamps';
            multi.durations{idx} = visuals.durations';

            multi.orth{idx} = 0; % do not orthogonalise them

            pix = 0;
            for i = 1:numel(fields)
                if ~ismember(fields{i}, {'timestamps', 'durations'}) 
                    if all(visuals.(fields{i}) == visuals.(fields{i})(1))
                        % constant
                        continue
                    end
                    pix = pix + 1;
                    multi.pmod(idx).name{pix} = fields{i};
                    multi.pmod(idx).param{pix} = visuals.(fields{i});
                    multi.pmod(idx).poly{pix} = 1;
                end
            end

            % GLM 8: on/off nuisance regressors
            %
            [onoff] = get_onoff(subj_id, run, conn, true);
            fields = fieldnames(onoff);
            for i = 1:numel(fields)
                idx = idx + 1;
                multi.names{idx} = fields{i};
                multi.onsets{idx} = onoff.(fields{i});
                multi.durations{idx} = zeros(size(multi.onsets{idx}));
            end


        % frame nuisance regressors: visual control regressors for each frame 
        % same as GLM 7 excep 1 per GLM b/c they're highly correlated and we get nothing
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        % TODO tightly coupled to get_visuals
        %
        case {10,11,12,13,14,15,16,17,18,19,20}  

            [fields, visuals] = get_visuals(subj_id, run, conn, true);
            fields

            multi.names{1} = 'frames';
            multi.onsets{1} = visuals.timestamps';
            multi.durations{1} = visuals.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            glm_idx = 9;
            for i = 1:numel(fields)
                if ~ismember(fields{i}, {'timestamps', 'durations'}) 
                    if all(visuals.(fields{i}) == visuals.(fields{i})(1))
                        % constant
                        continue
                    end

                    glm_idx = glm_idx + 1;
                    if glm_idx == glmodel
                        multi.pmod(1).name{1} = fields{i};
                        multi.pmod(1).param{1} = visuals.(fields{i});
                        multi.pmod(1).poly{1} = 1;
                    end
                end
            end

        % theory_change_flag + control regressors
        % i.e. GLM 3 + GLM 9 = 3, 1, 5, 7, 8
        % TODO dedupe
        %
        case 21

            idx = 0;

            % from GLM 3: theory_change_flag
            %
            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.theory_change_flag_onsets;


            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;

            % from GLM 1: game instance boxcar regressor
            %
            [game_names, onsets, durs] = get_games(subj_id, run, conn);
            for i = 1:numel(game_names)
                idx = idx + 1;
                multi.names{idx} = game_names{i};
                multi.onsets{idx} = onsets{i};
                multi.durations{idx} = durs{i};
            end

            % from GLM 5: keyholds nuisance regressors
            %
            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn, true);
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

            % GLM 7: frame nuisance regressors
            %
            [fields, visuals] = get_visuals(subj_id, run, conn, true);
            idx = idx + 1;
            multi.names{idx} = 'frames';
            multi.onsets{idx} = visuals.timestamps';
            multi.durations{idx} = visuals.durations';

            multi.orth{idx} = 0; % do not orthogonalise them

            pix = 0;
            for i = 1:numel(fields)
                if ~ismember(fields{i}, {'timestamps', 'durations'}) 
                    if all(visuals.(fields{i}) == visuals.(fields{i})(1))
                        % constant
                        continue
                    end
                    pix = pix + 1;
                    multi.pmod(idx).name{pix} = fields{i};
                    multi.pmod(idx).param{pix} = visuals.(fields{i});
                    multi.pmod(idx).poly{pix} = 1;
                end
            end

            % GLM 8: on/off nuisance regressors
            %
            [onoff] = get_onoff(subj_id, run, conn, true);
            fields = fieldnames(onoff);
            for i = 1:numel(fields)
                idx = idx + 1;
                multi.names{idx} = fields{i};
                multi.onsets{idx} = onoff.(fields{i});
                multi.durations{idx} = zeros(size(multi.onsets{idx}));
            end
    
        % momchil TODO FIXME: this is a bad idea; events are too close => must use actual deconvolution
        % beta series -- separate regressor for each TR for each level for each game
        % for RSA, effectively deconvolution
        % similar to GLM 23 in Exploration
        %
        case 22
            assert(false);

            idx = 0;

            EXPT = vgdl_expt();

            [game_names, onsets, durs] = get_games(subj_id, run, conn);

            for i = 1:numel(game_names)
                for j = 1:length(onsets{i})
                    ons = [onsets{i}(j):EXPT.TR:onsets{i}(j)+durs{i}(j)];
                    for k = 1:length(ons)
                        idx = idx + 1;
                        % impulse regressors tiling every instance at 2 s intervals 
                        multi.names{idx} = sprintf('run_%d_block_%d_instance_%d_time_%d', run_id, i, j, k);
                        multi.onsets{idx} = ons(k);
                        multi.durations{idx} = 0;
                    end
                end
            end

        % functional connectivity GLM 
        % TODO FIXME this doesn't work b/c we can't use beta series here, TRs are too close
        case 23
            assert(false);

            onsets = [];
            beta_series_glm = 22;

            % cheat: just get the onsets from the beta series GLM, b/c we can't do it otherwise on the cluster 
            % (and we can't pre-run it and cache the multi locally either, b/c the betas...)
            %
            multi_series = vgdl_create_multi(beta_series_glm, subj_id, run_id);
            for i = 1:length(multi_series.onsets)
                onsets = [onsets multi_series.onsets{i}];
            end

            EXPT = vgdl_expt();

            multi.names{1} = 'run';
            multi.onsets{1} = onsets;
            multi.durations{1} = zeros(size(onsets));

            R_IFG = ccnl_get_beta_series(EXPT, beta_series_glm, subj_id, sprintf('run_%d_', run_id), 'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii');
            R_IFG = mean(R_IFG,2);

            multi.pmod(1).name{1} = 'R_IFG';
            multi.pmod(1).param{1} = R_IFG';
            multi.pmod(1).poly{1} = 1;

        % 10 s boxcars tiling every instance, i.e. coarse beta series
        % for RSA, middle ground between GLM 1 (too coarse) and GLM 22 (too fine grained; doesn't fit)
        case 24

            idx = 0;

            [game_names, onsets, durs] = get_games(subj_id, run, conn);

            for i = 1:numel(game_names)
                for j = 1:length(onsets{i})
                    dt = (durs{i}(j) + 0.0001) / 6; % assumes instances are about 60 s long
                    ons = [onsets{i}(j) : dt : onsets{i}(j)+durs{i}(j)];
                    for k = 1:length(ons)
                        idx = idx + 1;
                        % 10 s boxcar regressors tiling every instance 
                        multi.names{idx} = sprintf('run_%d_block_%d_instance_%d_boxcar_%d', run_id, i, j, k);
                        multi.onsets{idx} = ons(k);
                        multi.durations{idx} = dt;
                    end
                end
            end


        % condition = game, boxcars over instances
        % for RSA, middle ground between GLM 1 and GLM 24
        %
        case 25

            idx = 0;

            [game_names, onsets, durs] = get_games(subj_id, run, conn);

            for i = 1:numel(game_names)
                % instance boxcar regressor
                for j = 1:length(onsets{i})
                    idx = idx + 1;
                    multi.names{idx} = sprintf('%s_run_%d_block_%d_instance_%d', game_names{i}, run_id, i, j);
                    multi.onsets{idx} = onsets{i}(j);
                    multi.durations{idx} = durs{i}(j);
                end
            end


        % individual theory-based regressors
        % same idea as GLM 10-20
        % TODO tight coupling with get_regressors.m
        %
        case {26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44}

            [regs, ~, fields] = get_regressors(subj_id, run, conn, true);
            fields

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            glm_idx = 25;
            for i = 1:numel(fields)
                if ~ismember(fields{i}, {'timestamps', 'durations'}) && ~contains(fields{i}, '_onsets')
                    glm_idx = glm_idx + 1;
                    if glm_idx == glmodel
                        idx = 0;

                        for j = 1:size(regs.(fields{i}), 2)
                            if all(regs.(fields{i})(:,j) == regs.(fields{i})(1,j))
                                % constant
                                continue
                            end

                            idx = idx + 1;
                            if size(regs.(fields{i}), 2) == 1
                                multi.pmod(1).name{idx} = fields{i};
                            else
                                multi.pmod(1).name{idx} = [fields{i}, '_', num2str(j)];
                            end
                            multi.pmod(1).param{idx} = regs.(fields{i})(:,j);
                            multi.pmod(1).poly{idx} = 1;
                        end

                    end
                end
            end


        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

    if save_output
        save(filename, 'multi', '-v7.3'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
    end

    if exist('conn', 'var')
        close(conn);
    end

