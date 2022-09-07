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
    %w

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
    run_id = get_behavioral_run_id(subj_id, SPM_run_id)
    fprintf('behavioral run_id %d \n', run_id);

    filename = fullfile(get_mat_dir(false), sprintf('vgdl_create_multi_glm%d_subj%d_run%d.mat', glmodel, subj_id, run_id));
    filename

    % hack to make it work on the cluster until they install MongoDB toolbox
    % just pre-generate the multi's locally, then load them on the cluster
    if ~ismember(glmodel, 23)
        try
            conn = mongo('holy7c22108.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54', 'UserName', 'reader', 'Password', 'parolatamadafaqa')
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

    multi.names = {};
    multi.onsets = {};
    multi.durations = {};

    % GLMs
    %
    switch glmodel

        % condition = game, boxcars over blocks
        % look for systematic differences across games for things like agency, etc
        % e.g. contrast: 'chase + helper - butterflies - aliens'
        %
        case 1 

            multi = add_games_to_multi(multi, subj_id, run, conn);

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

            regs = get_regressors(subj_id, run, conn, true);
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

            multi = add_keyholds_to_multi(multi, subj_id, run, conn);


        % key presses, (copied) subset of GLM 2, compare to GLM 5
        % nuisance regressors
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        %
        case 6

            multi = add_keypresses_to_multi(multi, subj_id, run, conn);


        % frame nuisance regressors: visual control regressors for each frame
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        %
        case 7

            multi = add_visuals_to_multi(multi, subj_id, run, conn);

        % on/off nuisance regressors: visual control regressors for start/end block, instance, play
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        %
        case 8

            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % condition = game, boxcars over blocks (GLM 1)
        % + nuisance regressors (buttons, frames, on/off) (GLMs 5, 7, 8)
        % look for systematic differences across games for things like agency, etc
        % e.g. contrast: 'chase + helper - butterflies - aliens'
        %
        case 9

            % from GLM 1: game instance boxcar regressor
            %
            multi = add_games_to_multi(multi, subj_id, run, conn);

            % from GLM 5: keyholds nuisance regressors
            %
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);

            % GLM 7: frame nuisance regressors
            %
            multi = add_visuals_to_multi(multi, subj_id, run, conn);

            % GLM 8: on/off nuisance regressors
            %
            multi = add_onoff_to_multi(multi, subj_id, run, conn);


        % frame nuisance regressors: visual control regressors for each frame 
        % same as GLM 7 excep 1 per GLM b/c they're highly correlated and we get nothing
        % this is exploratory, so ok to look at them separately and incorporate into other GLMs later
        % TODO tightly coupled to get_visuals
        %
        case {10,11,12,13,14,15,16,17,18,19,20}  

            [~, visuals] = get_visuals(subj_id, run, conn, true);
            fields = {
				'new_sprites', ...
				'killed_sprites', ...
				'sprites', ...
				'non_walls', ...
				'avatar_moved', ...
				'moved', ...
				'movable', ...
				'collisions', ...
				'effects', ...
				'sprite_groups', ...
				'changed'}
            field = fields{glmodel - 9};
            field

            multi.names{1} = 'frames';
            multi.onsets{1} = visuals.timestamps';
            multi.durations{1} = visuals.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            multi.pmod(1).name{1} = field;
            multi.pmod(1).param{1} = visuals.(field);
            multi.pmod(1).poly{1} = 1;

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

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

    
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


        % individual theory-based regressors, as pmods throughout the episodes
        % same idea as GLM 10-20
        % TODO tight coupling with get_regressors.m
        %
        case {26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,57,72,73,74}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true);
            fields = {
                'theory_change_flag', ...
                'sprite_change_flag', ...
                'termination_change_flag', ...
                'newEffects_flag', ...
                'likelihood', ...
                'sum_lik_play', ...
                'n_ts', ...
                'num_effects', ... # len(effectListByColor)
                'R_GG', ...
                'R_GGs', ...
                'R_SG', ...
                'R_SGs', ...
                'interaction_change_flag', ... % 38
                'S_len', ...
                'I_len', ...
                'T_len', ...
                'Igen_len', ... % 42
                'Tnov_len', ... % 43
                'Ip_len', ... % 44
                'dS_len', ... % 45
                'dI_len', ... % 46
                'dT_len', ... % 47
                'dIgen_len', ... % 48
                'dTnov_len', ... $ 49
                'dIp_len'}
            map = containers.Map(26:50, fields);
            map(57) = 'surprise';
            map(72) = 'replan_flag';
            map(73) = 'sum_lik';
            map(74) = 'spriteKL';

            field = map(glmodel);
            field

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            idx = 0;

            for j = 1:size(regs.(field), 2)
                if all(regs.(field)(:,j) == regs.(field)(1,j))
                    % constant
                    continue
                end

                idx = idx + 1;
                if size(regs.(field), 2) == 1
                    multi.pmod(1).name{idx} = field;
                else
                    multi.pmod(1).name{idx} = [field, '_', num2str(j)];
                end
                multi.pmod(1).param{idx} = regs.(field)(:,j);
                multi.pmod(1).poly{idx} = 1;

                multi.pmod(1).param{idx}(isnan(multi.pmod(1).param{idx})) = 0; % TODO happens for lik during first few timesteps; ideally, remove altogether
            end

        % interaction_change_flag 
        % c/p GLM 3
        %
        case 51

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.interaction_change_flag_onsets;

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


        % termination_change_flag 
        % c/p GLM 3
        %
        case 52

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.termination_change_flag_onsets;

            if length(regs.termination_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

        % sprite_change_flag, interaction_change_flag, termination_change_flag 
        % union GLMs 51..52 + sprites
        %
        case 53

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % sometimes interactions and terminations are identical
            if length(regs.termination_change_flag_onsets) > 0 && (length(regs.termination_change_flag_onsets) ~= length(regs.interaction_change_flag_onsets) || any(regs.interaction_change_flag_onsets ~= regs.termination_change_flag_onsets))
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


        % sprite_change_flag, interaction_change_flag, termination_change_flag as pmods
        % GLM 53 but pmods
        %
        case 54

            idx = 0;

            [regs, ~, fields] = get_regressors(subj_id, run, conn, true);
            fields

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            if any(regs.sprite_change_flag)
                idx = idx + 1;
                multi.pmod(1).name{idx} = 'sprite_change_flag';
                multi.pmod(1).param{idx} = regs.sprite_change_flag;
                multi.pmod(1).poly{idx} = 1;
            end

            if any(regs.interaction_change_flag)
                idx = idx + 1;
                multi.pmod(1).name{idx} = 'interaction_change_flag';
                multi.pmod(1).param{idx} = regs.interaction_change_flag;
                multi.pmod(1).poly{idx} = 1;
            end

            % sometimes interactions and terminations are identical
            if any(regs.termination_change_flag) && any(xor(regs.interaction_change_flag, regs.termination_change_flag))
                idx = idx + 1;
                multi.pmod(1).name{idx} = 'termination_change_flag';
                multi.pmod(1).param{idx} = regs.termination_change_flag;
                multi.pmod(1).poly{idx} = 1;
            end


        % theory_change_flag - avatar_collision_flag (orth)
        % #Sam
        %
        case 55

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            tc = regs.theory_change_flag;
            tc = logical(tc);
            assert(all(regs.timestamps(tc) == regs.theory_change_flag_onsets));

            [~, visuals] = get_visuals(subj_id, run, conn, true);
            av = visuals.avatar_collision_flag;
            av = logical(av);
            which = ismember(visuals.timestamps, regs.state_timestamps); % cross-reference them; annoying stuff
            av = av(which);

            av(tc) = 0; % orthogonalise w.r.t. theory_change_flag

            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = regs.timestamps(tc);
            multi.durations{idx} = zeros(size(multi.onsets{idx}));

            idx = idx + 1;
            multi.names{idx} = 'avatar_collision_flag';
            multi.onsets{idx} = regs.timestamps(av);
            multi.durations{idx} = zeros(size(multi.onsets{idx}));


        % theory_change_flag - collisions (orth)
        % similar to 55
        %
        case 56

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            tc = regs.theory_change_flag;
            tc = logical(tc);
            assert(all(regs.timestamps(tc) == regs.theory_change_flag_onsets));

            [~, visuals] = get_visuals(subj_id, run, conn, true);
            ef = visuals.effectsByCol; % same as newTimeStep_flag 
            ef = ef > 0;
            which = ismember(visuals.timestamps, regs.state_timestamps); % cross-reference them; annoying stuff
            ef = ef(which);
            %assert(all(ef == regs.newTimeStep_flag)); fails TODO weird

            ef(tc) = 0; % orthogonalize w.r.t. theory_change_flag

            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = regs.timestamps(tc);
            multi.durations{idx} = zeros(size(multi.onsets{idx}));

            idx = idx + 1;
            multi.names{idx} = 'collision_flag';
            multi.onsets{idx} = regs.timestamps(ef);
            multi.durations{idx} = zeros(size(multi.onsets{idx}));

        % case 57 is surprise
        % theory_change_flag vs. surprise 
        %
        case 58

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            tc = regs.theory_change_flag;
            tc = logical(tc);
            assert(all(regs.timestamps(tc) == regs.theory_change_flag_onsets));

            su = regs.surprise;
            su(isnan(su)) = 0;
            su = logical(su);

            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = regs.timestamps(tc);
            multi.durations{idx} = zeros(size(multi.onsets{idx}));

            idx = idx + 1;
            multi.names{idx} = 'surprise_flag';
            multi.onsets{idx} = regs.timestamps(su);
            multi.durations{idx} = zeros(size(multi.onsets{idx}));


        % old theory_change_flag 
        % same as GLM 3 but with regressors_2020_04_01_fullTSList_nonrefactor
        %
        case 59

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true, 'regressors_1');
            onsets = regs.theory_change_flag_onsets;

            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;

        % old theory_change_flag 
        % same as GLM 3 but with regressors_2020_05_19_finalTimeStep_reset_nonrefactored
        %
        case 60

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true, 'regressors_2');
            onsets = regs.theory_change_flag_onsets;

            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;

        % old theory_change_flag 
        % same as GLM 3 but with regressors_2020_05_24_finalTimeStep_reset_refactored
        %
        case 61

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true, 'regressors_3');
            onsets = regs.theory_change_flag_onsets;

            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;


        % old theory_change_flag + control regressors
        % same as GLM 21 but with regressors_2020_04_01_fullTSList_nonrefactor
        %
        case 62

            idx = 0;

            % from GLM 3: theory_change_flag
            %
            regs = get_regressors(subj_id, run, conn, true, 'regressors_1');
            onsets = regs.theory_change_flag_onsets;


            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);


        % old theory_change_flag + control regressors
        % same as GLM 21 but with regressors_2020_05_19_finalTimeStep_reset_nonrefactored
        %
        case 63

            idx = 0;

            % from GLM 3: theory_change_flag
            %
            regs = get_regressors(subj_id, run, conn, true, 'regressors_2');
            onsets = regs.theory_change_flag_onsets;


            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);


        % old theory_change_flag + control regressors
        % same as GLM 21 but with regressors_2020_05_24_finalTimeStep_reset_refactored
        %
        case 64

            idx = 0;

            % from GLM 3: theory_change_flag
            %
            regs = get_regressors(subj_id, run, conn, true, 'regressors_3');
            onsets = regs.theory_change_flag_onsets;


            idx = idx + 1;
            multi.names{idx} = 'theory_change_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;


            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % old sprite_change_flag, interaction_change_flag, termination_change_flag 
        % same as GLM 53 but with regressors_2020_04_01_fullTSList_nonrefactor

        %
        case 65

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true, 'regressors_1');

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            if length(regs.termination_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


        % old sprite_change_flag, interaction_change_flag, termination_change_flag 
        % same as GLM 53 but with regressors_2020_05_19_finalTimeStep_reset_nonrefactored

        %
        case 66

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true, 'regressors_2');

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.termination_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


        % interaction_change_flag or termination_change_flag 
        % similar idea to GLM 3 except without the sprite changes 
        % also looks like sometimes theory_change_flag = 1 but nothing actually chanfed, just some interactions got shuffled around
        %
        case 67

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            ic = regs.interaction_change_flag;
            ic = logical(ic);
            assert(all(regs.timestamps(ic) == regs.interaction_change_flag_onsets));

            tec = regs.termination_change_flag;
            tec = logical(tec);
            assert(all(regs.timestamps(tec) == regs.termination_change_flag_onsets));

            if any(ic | tec)
                idx = idx + 1;
                multi.names{idx} = 'interaction_or_termination_change_flag';
                multi.onsets{idx} = regs.timestamps(ic | tec);
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            %{
        % old sprite_change_flag, interaction_change_flag, termination_change_flag 
        % same as 53 but with regressors_2020_05_24_finalTimeStep_reset_refactored
        %
        case 68

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true, 'regressors_3');

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            idx = idx + 1;
            multi.names{idx} = 'interaction_change_flag';
            multi.onsets{idx} = regs.interaction_change_flag_onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;


            idx = idx + 1;
            multi.names{idx} = 'termination_change_flag';
            multi.onsets{idx} = regs.termination_change_flag_onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;


        % old sprite_change_flag, interaction_change_flag, termination_change_flag as pmods
        % same as GLM 54 but with regressors_2020_05_24_finalTimeStep_reset_refactored
        %
        case 69

            idx = 0;

            [regs, ~, fields] = get_regressors(subj_id, run, conn, true, 'regressors_3');
            fields

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            if any(regs.sprite_change_flag)
                idx = idx + 1;
                multi.pmod(1).name{idx} = 'sprite_change_flag';
                multi.pmod(1).param{idx} = regs.sprite_change_flag;
                multi.pmod(1).poly{idx} = 1;
            end

            idx = idx + 1;
            multi.pmod(1).name{idx} = 'interaction_change_flag';
            multi.pmod(1).param{idx} = regs.interaction_change_flag;
            multi.pmod(1).poly{idx} = 1;

            idx = idx + 1;
            multi.pmod(1).name{idx} = 'termination_change_flag';
            multi.pmod(1).param{idx} = regs.termination_change_flag;
            multi.pmod(1).poly{idx} = 1;

            %}


        % sprite_change_flag, interaction_change_flag, termination_change_flag 
        %  + control regressors
        % merge of GLM 21 and GLM 53
        % i.e. GLM 53 + GLM 9 = 53, 1, 5, 7, 8
        %
        case 68

            idx = 0;

            % from GLM 53: sprite_change_flag, interaction_change_flag, termination_change_flag 

            %
            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            if length(regs.termination_change_flag_onsets) > 0 && (length(regs.termination_change_flag_onsets) ~= length(regs.interaction_change_flag_onsets) || any(regs.interaction_change_flag_onsets ~= regs.termination_change_flag_onsets))
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % surprise, likelihood, sum_lik_play
        %  + control regressors
        % see GLM 21 
        %
        case 69

            % from GLM 30,31,32
            fields = {
                'likelihood', ...
                'sum_lik_play', ...
                'surprise'}
            multi = add_frames_with_pmods(struct, subj_id, run, conn, fields)

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % sum_lik_play vs. n_ts
        %
        case 70

            % from GLM 30,31,32,57: sprite_change_flag, interaction_change_flag, termination_change_flag 

            %
            fields = {
                'sum_lik_play', ...
                'n_ts' ...
                }
            multi = add_frames_with_pmods(struct, subj_id, run, conn, fields)

        % replan_flag 
        %
        case 71

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.replan_flag_onsets;

            idx = idx + 1;
            multi.names{idx} = 'replan_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;

        % interaction_change_flag or termination_change_flag, with nuisance regressors
        % i.e. GLM 67 + GLM 9
        % similar idea to GLM 21
        %
        case 75

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            ic = regs.interaction_change_flag;
            ic = logical(ic);
            assert(all(regs.timestamps(ic) == regs.interaction_change_flag_onsets));

            tec = regs.termination_change_flag;
            tec = logical(tec);
            assert(all(regs.timestamps(tec) == regs.termination_change_flag_onsets));

            if any(ic | tec)
                idx = idx + 1;
                multi.names{idx} = 'interaction_or_termination_change_flag';
                multi.onsets{idx} = regs.timestamps(ic | tec);
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % spriteKL
        %  + control regressors
        % like GLM 21; also see GLM 69
        %
        case 76

            % from GLM 30,31,32,57: sprite_change_flag, interaction_change_flag, termination_change_flag 

            %
            fields = {
                'spriteKL', ...
                }
            multi = add_frames_with_pmods(struct, subj_id, run, conn, fields)

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);


        % GLM 9 but without game id
        % i.e. 5, 7, 8
        %
        case 77

            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);
   

        % condition = game, boxcars over blocks (GLM 1)
        % + key presses (GLM 5)
        % 
        case 78

            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);

        % surprise_flag 
        % see 58
        %
        case 79

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            su = regs.surprise;
            su(isnan(su)) = 0;
            su = logical(su);
            onsets = regs.timestamps(su);

            idx = idx + 1;
            multi.names{idx} = 'surprise_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));


        % surprise_flag + control regressors
        % like GLM 21
        % i.e. GLM 79 + GLM 9 = 79, 1, 5, 7, 8
        %
        case 80

            idx = 0;

            % from GLM 79: surprise_flag 
            %
            regs = get_regressors(subj_id, run, conn, true);

            su = regs.surprise;
            su(isnan(su)) = 0;
            su = logical(su);
            onsets = regs.timestamps(su);

            idx = idx + 1;
            multi.names{idx} = 'surprise_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);


        % replan_flag + control regressors
        % like GLM 21
        % i.e. GLM 71 + GLM 9 = 79, 1, 5, 7, 8
        %
        case 81

            idx = 0;

            % from GLM 79: surprise_flag 
            %
            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.replan_flag_onsets;

            idx = idx + 1;
            multi.names{idx} = 'replan_flag';
            multi.onsets{idx} = onsets;
            multi.durations{idx} = zeros(size(multi.onsets{idx}));;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

    
        % interaction_change_flag 
        %  + control regressors
        % like GLM 21
        % i.e. GLM 53 + GLM 9 = 53, 1, 5, 7, 8
        %
        case 82

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % termination_change_flag 
        %  + control regressors
        % like GLM 21
        % i.e. GLM 53 + GLM 9 = 53, 1, 5, 7, 8
        %
        case 83

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.termination_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % likelihood
        %  + control regressors
        % see GLM 21 and/or 69
        %
        case 84

            multi = add_frames_with_pmods(struct, subj_id, run, conn, {'likelihood'})

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % sprite_change_flag
        % c/p GLM 3
        %
        case 85

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.sprite_change_flag_onsets;

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


        % sprite_change_flag + nuisance regressors
        case 86

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % interaction_change_flag, termination_change_flag 
        % union GLMs 51..52
        %
        case 87

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % sometimes interactions and terminations are identical
            if length(regs.termination_change_flag_onsets) > 0 && (length(regs.termination_change_flag_onsets) ~= length(regs.interaction_change_flag_onsets) || any(regs.interaction_change_flag_onsets ~= regs.termination_change_flag_onsets))
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

        % interaction_change_flag, termination_change_flag 
        % union GLMs 51..52 + nuisance regressors
        %
        case 88

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % sometimes interactions and terminations are identical
            if length(regs.termination_change_flag_onsets) > 0 && (length(regs.termination_change_flag_onsets) ~= length(regs.interaction_change_flag_onsets) || any(regs.interaction_change_flag_onsets ~= regs.termination_change_flag_onsets))
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % spriteKL_flag, i.e. discretized (so no pmod; b/c sprite_change_flag is a little flaky)
        % TODO doesn't really work, some runs end up empty (e.g. subj 23, run 3)
        % see also GLM 75
        %
        case 89

            idx = 0;
            multi = struct;

            regs = get_regressors(subj_id, run, conn, true);
            sc = logical(regs.spriteKL > 0);

            %if any(sc)
                idx = idx + 1;
                multi.names{idx} = 'spriteKL_flag'
                multi.onsets{idx} = regs.timestamps(sc);
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            %end

        % spriteKL_flag + nuisance regressors        % see also GLM 75
        % TODO confounded with new_sprites
        %
        case 90

            idx = 0;
            multi = struct;

            regs = get_regressors(subj_id, run, conn, true);
            sc = logical(regs.spriteKL > 0);

            if any(sc)
                idx = idx + 1;
                multi.names{idx} = 'spriteKL_flag'
                multi.onsets{idx} = regs.timestamps(sc);
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % sprite_change_flag or interaction_change_flag
        % see GLM 67
        %
        case 91

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            sc = regs.sprite_change_flag;
            sc = logical(sc);
            assert(all(regs.timestamps(sc) == regs.sprite_change_flag_onsets));

            ic = regs.interaction_change_flag;
            ic = logical(ic);
            assert(all(regs.timestamps(ic) == regs.interaction_change_flag_onsets));

            if any(sc | ic)
                idx = idx + 1;
                multi.names{idx} = 'sprite_or_interaction_change_flag';
                multi.onsets{idx} = regs.timestamps(sc | ic);
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

        % sprite_change_flag or termination_change_flag
        % see GLM 67
        %
        case 92

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            sc = regs.sprite_change_flag;
            sc = logical(sc);
            assert(all(regs.timestamps(sc) == regs.sprite_change_flag_onsets));

            tec = regs.termination_change_flag;
            tec = logical(tec);
            assert(all(regs.timestamps(tec) == regs.termination_change_flag_onsets));

            if any(sc | tec)
                idx = idx + 1;
                multi.names{idx} = 'sprite_or_termination_change_flag';
                multi.onsets{idx} = regs.timestamps(sc | tec);
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

        % sprite_change_flag or interaction_change_flag
        % + nuisance regressors
        % i.e. GLM 91 + GLM 9
        %
        case 93

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            sc = regs.sprite_change_flag;
            sc = logical(sc);
            assert(all(regs.timestamps(sc) == regs.sprite_change_flag_onsets));

            ic = regs.interaction_change_flag;
            ic = logical(ic);
            assert(all(regs.timestamps(ic) == regs.interaction_change_flag_onsets));

            if any(sc | ic)
                idx = idx + 1;
                multi.names{idx} = 'sprite_or_interaction_change_flag';
                multi.onsets{idx} = regs.timestamps(sc | ic);
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % sprite_change_flag or termination_change_flag
        % + nuisance regressors
        % i.e. GLM 92 + GLM 9
        %
        case 94

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);
            sc = regs.sprite_change_flag;
            sc = logical(sc);
            assert(all(regs.timestamps(sc) == regs.sprite_change_flag_onsets));

            tec = regs.termination_change_flag;
            tec = logical(tec);
            assert(all(regs.timestamps(tec) == regs.termination_change_flag_onsets));

            if any(sc | tec)
                idx = idx + 1;
                multi.names{idx} = 'sprite_or_termination_change_flag';
                multi.onsets{idx} = regs.timestamps(sc | tec);
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % sprite_change_flag, interaction_change_flag 
        % see GLM 53
        %
        case 95

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

        % sprite_change_flag, termination_change_flag 
        % see GLM 53
        %
        case 96

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.termination_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

        % sprite_change_flag, interaction_change_flag 
        % + nuisance regressors
        % see GLM 68
        %
        case 97

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % sprite_change_flag, termination_change_flag 
        % + nuisance regressors
        % see GLM 68
        %
        case 98

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.termination_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % # exploratory goals: generic interactions and novelty rule terminations 
        % see GLMs 42 and 43
        %
        case 99

            fields = {
                'Igen_len', ...
                'Tnov_len'}
            multi = add_frames_with_pmods(struct, subj_id, run, conn, fields)

        % # exploratory goals: generic interactions and novelty rule terminations 
        % + nuisance regressors
        % see GLM 99
        %
        case 100

            fields = {
                'Igen_len', ...
                'Tnov_len'}
            multi = add_frames_with_pmods(struct, subj_id, run, conn, fields)

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn);

        % like GLM 9, but with play on/off
        %
        case 101

            % from GLM 1: game instance boxcar regressor
            %
            multi = add_games_to_multi(multi, subj_id, run, conn);

            % from GLM 5: keyholds nuisance regressors
            %
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);

            % GLM 7: frame nuisance regressors
            %
            multi = add_visuals_to_multi(multi, subj_id, run, conn);

            % GLM 8: on/off nuisance regressors
            %
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % like GLM 21, but with play on/off
        %
        case 102

            % from GLM 3: theory_change_flag
            %
            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.theory_change_flag_onsets;

            multi.names{1} = 'theory_change_flag';
            multi.onsets{1} = onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 86, but with play on/off
        case 103

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 82, but with play on/off
        case 104

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 83, but with play on/off
        case 105

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.termination_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 68, but with play on/off
        case 106

            idx = 0;

            % from GLM 53: sprite_change_flag, interaction_change_flag, termination_change_flag 

            %
            regs = get_regressors(subj_id, run, conn, true);

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            if length(regs.termination_change_flag_onsets) > 0 && (length(regs.termination_change_flag_onsets) ~= length(regs.interaction_change_flag_onsets) || any(regs.interaction_change_flag_onsets ~= regs.termination_change_flag_onsets))
                idx = idx + 1;
                multi.names{idx} = 'termination_change_flag';
                multi.onsets{idx} = regs.termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % .State features
        case 107

            fields = {'sprites', 'new_sprites', 'killed_sprites', 'collisions', 'effects', 'sprite_groups', 'effectsByCol', 'win', 'loss', 'ended', 'score', 'dscore'};
            multi = add_visuals_to_multi(multi, subj_id, run, conn, fields);

        % .irrelevant features
        case 108

            fields = {'non_walls', 'avatar_moved', 'moved', 'movable', 'changed', 'avatar_collision_flag'};
            multi = add_visuals_to_multi(multi, subj_id, run, conn, fields);

        % Beta series GLM for theory_change_flag
        case {109, 120, 121, 122}

            minimum_onset_separation = 10; % s

            % optionally offset from the theory_change_flag time
            switch glmodel
                case 109
                    lag = 0;
                case 120
                    lag = 2; % 2 s = 1 TR
                case 121
                    lag = 4; % 4 s = 2 TRs
                case 122
                    lag = 6; % 6 s = 3 TRs
                otherwise
                    assert(false);
            end

            % from GLM 3: theory_change_flag
            %
            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.theory_change_flag_onsets;

            idx = 0;
            last_onset = -1000;
            for i = 1:length(onsets)
                if onsets(i) - last_onset >= minimum_onset_separation
                    idx = idx + 1;
                    multi.names{idx} = sprintf('theory_change_flag_%d', i);
                    multi.onsets{idx} = onsets(i) + lag;
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));
                    last_onset = onsets(i);
                end
            end

        % replan_flag + keypresses
        %
        case 110

            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.replan_flag_onsets;

            multi.names{1} = 'replan_flag';
            multi.onsets{1} = onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            % nuisance regressors
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);

        % replan_flag + keypresses + on/off
        %
        case 111

            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.replan_flag_onsets;

            multi.names{1} = 'replan_flag';
            multi.onsets{1} = onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            % nuisance regressors
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % number of exploratory goals
        case 112

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true);
            exploratory_goals = regs.Igen_len + regs.Tnov_len;

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            multi.pmod(1).name{1} = 'exploratory_goals';
            multi.pmod(1).param{1} = exploratory_goals;
            multi.pmod(1).poly{1} = 1;

        % number of exploitative goals
        case 113

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true);
            exploitative_goals = regs.T_len;

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            if std(exploitative_goals) > 1e-6
                % make sure they're not constant e.g. subj 24 run 3
                multi.pmod(1).name{1} = 'exploitative_goals';
                multi.pmod(1).param{1} = exploitative_goals;
                multi.pmod(1).poly{1} = 1;
            end

        % number of exploratory goals + keypresses + on/off
        case 114

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true);
            exploratory_goals = regs.Igen_len + regs.Tnov_len;

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            multi.pmod(1).name{1} = 'exploratory_goals';
            multi.pmod(1).param{1} = exploratory_goals;
            multi.pmod(1).poly{1} = 1;

            %  nuisance regressors
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % number of exploitative goals + keypresses + on/off
        case 115

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true);
            exploitative_goals = regs.T_len;

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            if std(exploitative_goals) > 1e-6
                multi.pmod(1).name{1} = 'exploitative_goals';
                multi.pmod(1).param{1} = exploitative_goals;
                multi.pmod(1).poly{1} = 1;
            end

            %  nuisance regressors
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % exploration versus exploitation, boxcars
        case {116, 117, 118, 119}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps); % cross-reference them; annoying stuff
            if ~all(which)
                % TODO this seems outdated, never hits; remove eventually
                keyboard
            end
            acf = acf(which);
            ks = ks(which);
            ps = ps(which);
            pe = pe(which);
            win = win(which);

            %acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            %acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (tcf & acf);

            % initiation point for action sequences, both exploration and exploitation
            switch glmodel
                case {116, 118}
                    starts = exploit | explore | ps;
                case {117, 119}
                    starts = exploit | explore | ps | acf; % TODO rm not very useful
                otherwise
                    assert(false);
            end

            % convert to timestamps
            exploit_ts = regs.timestamps(exploit);
            explore_ts = regs.timestamps(explore);
            starts_ts = regs.timestamps(starts);

            idx = 0;

            % add exploitation boxcars
            if any(exploit_ts)
                idx = idx + 1;
                multi.names{idx} = 'exploit';
                multi.onsets{idx} = nan(size(exploit_ts));
                multi.durations{idx} = nan(size(multi.onsets{idx}));

                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ts)
                    en = exploit_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx}(i) = st;
                    multi.durations{idx}(i) = en - st;
                end
            end

            % add exploration boxcars
            if any(explore_ts)
                idx = idx + 1;
                multi.names{idx} = 'explore';
                multi.onsets{idx} = nan(size(explore_ts));
                multi.durations{idx} = nan(size(multi.onsets{idx}));

                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ts)
                    en = explore_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx}(i) = st;
                    multi.durations{idx}(i) = en - st;
                end
            end

            % optional control regressors
            switch glmodel
                case {116, 117}
                    disp('no control');

                case {118, 119}
                    if any(regs.theory_change_flag_onsets)
                        idx = idx + 1;
                        multi.names{idx} = 'theory_change_flag';
                        multi.onsets{idx} = regs.theory_change_flag_onsets;
                        multi.durations{idx} = zeros(size(multi.onsets{idx}));;
                    end

                otherwise
                    assert(false);
            end


        % number of exploratory goals & exploitative goals
        case {123, 124}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true);
            exploratory_goals = regs.Igen_len + regs.Tnov_len;
            exploitative_goals = regs.T_len;

            multi.names{1} = 'frames';
            multi.onsets{1} = regs.timestamps';
            multi.durations{1} = regs.durations';

            multi.orth{1} = 0; % do not orthogonalise them

            multi.pmod(1).name{1} = 'exploratory_goals';
            multi.pmod(1).param{1} = exploratory_goals;
            multi.pmod(1).poly{1} = 1;

            if std(exploitative_goals) > 1e-6
                % make sure they're not constant e.g. subj 24 run 3
                multi.pmod(1).name{2} = 'exploitative_goals';
                multi.pmod(1).param{2} = exploitative_goals;
                multi.pmod(1).poly{2} = 1;
            end

            % optionally add nuisance regressors
            if glmodel == 124
                multi = add_keyholds_to_multi(multi, subj_id, run, conn);
                multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});
            end


        % exploration versus exploitation, time to subgoal
        case {125, 126, 127, 128}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps);
            assert(all(which)) % Make sure they are all matching; we used to have to cross-reference them, but not any longer; this is just sent a check to ensure this stays the case

            acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (tcf & acf);

            % initiation point for action sequences, both exploration and exploitation
            starts = exploit | explore | ps;

            % convert to timestamps
            exploit_ix = find(exploit);
            explore_ix = find(explore);
            starts_ix = find(starts);

            idx = 0;

            idx = idx + 1;
            multi.names{idx} = 'frames';
            multi.onsets{idx} = regs.timestamps';
            multi.durations{idx} = regs.durations';

            multi.orth{idx} = 0; % do not orthogonalise them

            pidx = 0;

            % add exploitation ramps 
            if any(exploit_ix)
                ttsg = zeros(size(exploit));
                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ix)
                    eix = exploit_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    ttsg(six:eix) = regs.timestamps(eix) - regs.timestamps(six:eix);
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'exploit';
                multi.pmod(idx).param{pidx} = ttsg;
                multi.pmod(idx).poly{pidx} = 1;
            end

            % add exploration ramps 
            if any(explore_ix)
                ttsg = zeros(size(explore));
                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ix)
                    eix = explore_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    ttsg(six:eix) = regs.timestamps(eix) - regs.timestamps(six:eix);
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'explore';
                multi.pmod(idx).param{pidx} = ttsg;
                multi.pmod(idx).poly{pidx} = 1;
            end


            % optional control regressors
            if ismember(glmodel, [126, 128])
                if any(regs.theory_change_flag_onsets)
                    idx = idx + 1;
                    multi.names{idx} = 'theory_change_flag';
                    multi.onsets{idx} = regs.theory_change_flag_onsets;
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));;
                end
            end
            if ismember(glmodel, [127, 128])
                multi = add_keypresses_to_multi(multi, subj_id, run, conn);
            end


        % exploration versus exploitation, ramps
        case {129, 130, 131, 132}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps);
            assert(all(which)) % Make sure they are all matching; we used to have to cross-reference them, but not any longer; this is just sent a check to ensure this stays the case

            %acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            %acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (tcf & acf);

            % initiation point for action sequences, both exploration and exploitation
            starts = exploit | explore | ps;

            % convert to timestamps
            exploit_ix = find(exploit);
            explore_ix = find(explore);
            starts_ix = find(starts);

            idx = 0;

            idx = idx + 1;
            multi.names{idx} = 'frames';
            multi.onsets{idx} = regs.timestamps';
            multi.durations{idx} = regs.durations';

            multi.orth{idx} = 0; % do not orthogonalise them

            pidx = 0;

            % add exploitation ramps 
            if any(exploit_ix)
                ramp = zeros(size(exploit));
                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ix)
                    eix = exploit_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    ramp(six:eix) = linspace(0, 1, eix - six + 1);
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'exploit';
                multi.pmod(idx).param{pidx} = ramp;
                multi.pmod(idx).poly{pidx} = 1;
            end

            % add exploration ramps 
            if any(explore_ix)
                ramp = zeros(size(explore));
                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ix)
                    eix = explore_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    ramp(six:eix) = linspace(0, 1, eix - six + 1);
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'explore';
                multi.pmod(idx).param{pidx} = ramp;
                multi.pmod(idx).poly{pidx} = 1;
            end


            % optional control regressors
            if ismember(glmodel, [130, 132])
                if any(regs.theory_change_flag_onsets)
                    idx = idx + 1;
                    multi.names{idx} = 'theory_change_flag';
                    multi.onsets{idx} = regs.theory_change_flag_onsets;
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));;
                end
            end
            if ismember(glmodel, [131, 132])
                multi = add_keypresses_to_multi(multi, subj_id, run, conn);
            end


        % exploration versus exploitation, downward ramps
        case {133, 134, 135, 136}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps);
            assert(all(which)) % Make sure they are all matching; we used to have to cross-reference them, but not any longer; this is just sent a check to ensure this stays the case

            %acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            %acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (tcf & acf);

            % initiation point for action sequences, both exploration and exploitation
            starts = exploit | explore | ps;

            % convert to timestamps
            exploit_ix = find(exploit);
            explore_ix = find(explore);
            starts_ix = find(starts);

            idx = 0;

            idx = idx + 1;
            multi.names{idx} = 'frames';
            multi.onsets{idx} = regs.timestamps';
            multi.durations{idx} = regs.durations';

            multi.orth{idx} = 0; % do not orthogonalise them

            pidx = 0;

            % add exploitation ramps 
            if any(exploit_ix)
                ramp = zeros(size(exploit));
                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ix)
                    eix = exploit_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    ramp(six:eix) = linspace(1, 0, eix - six + 1);
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'exploit';
                multi.pmod(idx).param{pidx} = ramp;
                multi.pmod(idx).poly{pidx} = 1;
            end

            % add exploration ramps 
            if any(explore_ix)
                ramp = zeros(size(explore));
                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ix)
                    eix = explore_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    ramp(six:eix) = linspace(1, 0, eix - six + 1);
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'explore';
                multi.pmod(idx).param{pidx} = ramp;
                multi.pmod(idx).poly{pidx} = 1;
            end


            % optional control regressors
            if ismember(glmodel, [134, 136])
                if any(regs.theory_change_flag_onsets)
                    idx = idx + 1;
                    multi.names{idx} = 'theory_change_flag';
                    multi.onsets{idx} = regs.theory_change_flag_onsets;
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));;
                end
            end
            if ismember(glmodel, [135, 136])
                multi = add_keypresses_to_multi(multi, subj_id, run, conn);
            end


        % exploration versus exploitation, # keypresses to subgoal
        case {137, 138, 139, 140}

            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn, true);

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps);
            assert(all(which)) % Make sure they are all matching; we used to have to cross-reference them, but not any longer; this is just sent a check to ensure this stays the case

            %acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            %acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (tcf & acf);
            exploit = exploit & ~explore; % optionally orthogonalize them 

            % initiation point for action sequences, both exploration and exploitation
            starts = exploit | explore | ps;

            % convert to timestamps
            exploit_ix = find(exploit);
            explore_ix = find(explore);
            starts_ix = find(starts);

            idx = 0;

            idx = idx + 1;
            multi.names{idx} = 'frames';
            multi.onsets{idx} = regs.timestamps';
            multi.durations{idx} = regs.durations';

            multi.orth{idx} = 0; % do not orthogonalise them

            pidx = 0;

            % add exploitation ramps 
            if any(exploit_ix)
                ktsg = zeros(size(exploit));
                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ix)
                    eix = exploit_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    for j = six:eix
                        for k = 1:length(keypresses)
                            ktsg(j) = ktsg(j) + sum((keypresses{k} <= regs.timestamps(eix)) & (keypresses{k} >= regs.timestamps(j)));
                        end
                    end
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'exploit';
                multi.pmod(idx).param{pidx} = ktsg;
                multi.pmod(idx).poly{pidx} = 1;
            end

            % add exploration ramps 
            if any(explore_ix)
                ktsg = zeros(size(explore));
                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ix)
                    eix = explore_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    for j = six:eix
                        for k = 1:length(keypresses)
                            ktsg(j) = ktsg(j) + sum((keypresses{k} <= regs.timestamps(eix)) & (keypresses{k} >= regs.timestamps(j)));
                        end
                    end
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'explore';
                multi.pmod(idx).param{pidx} = ktsg;
                multi.pmod(idx).poly{pidx} = 1;
            end


            % optional control regressors
            if ismember(glmodel, [137, 140])
                if any(regs.theory_change_flag_onsets)
                    idx = idx + 1;
                    multi.names{idx} = 'theory_change_flag';
                    multi.onsets{idx} = regs.theory_change_flag_onsets;
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));;
                end
            end
            if ismember(glmodel, [139, 140])
                multi = add_keypresses_to_multi(multi, subj_id, run, conn);
            end


        % exploration versus exploitation, boxcars (see 116), modified
        case {141, 142, 143, 144}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            icf = logical(regs.interaction_change_flag);
            tecf = logical(regs.termination_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps); % cross-reference them; annoying stuff
            assert(all(which)) % Make sure they are all matching; we used to have to cross-reference them, but not any longer; this is just sent a check to ensure this stays the case

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (icf & acf) | (tecf & acf);
            exploit = exploit & ~explore; % orthogonalize them

            % initiation point for action sequences, both exploration and exploitation
            starts = exploit | explore | ps;

            % convert to timestamps
            exploit_ts = regs.timestamps(exploit);
            explore_ts = regs.timestamps(explore);
            starts_ts = regs.timestamps(starts);

            idx = 0;

            % add exploitation boxcars
            if any(exploit_ts)
                idx = idx + 1;
                multi.names{idx} = 'exploit';
                multi.onsets{idx} = nan(size(exploit_ts));
                multi.durations{idx} = nan(size(multi.onsets{idx}));

                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ts)
                    en = exploit_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx}(i) = st;
                    multi.durations{idx}(i) = en - st;
                end
            end

            % add exploration boxcars
            if any(explore_ts)
                idx = idx + 1;
                multi.names{idx} = 'explore';
                multi.onsets{idx} = nan(size(explore_ts));
                multi.durations{idx} = nan(size(multi.onsets{idx}));

                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ts)
                    en = explore_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx}(i) = st;
                    multi.durations{idx}(i) = en - st;
                end
            end

            % optional control regressors
            if ismember(glmodel, [142, 144])
                if any(regs.theory_change_flag_onsets)
                    idx = idx + 1;
                    multi.names{idx} = 'theory_change_flag';
                    multi.onsets{idx} = regs.theory_change_flag_onsets;
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));;
                end
            end
            if ismember(glmodel, [143, 144])
                multi = add_keypresses_to_multi(multi, subj_id, run, conn);
            end

        % exploration versus exploitation, downward ramps, orthogonalized
        case {145, 146, 147, 148}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            icf = logical(regs.interaction_change_flag);
            tecf = logical(regs.termination_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps);
            assert(all(which)) % Make sure they are all matching; we used to have to cross-reference them, but not any longer; this is just sent a check to ensure this stays the case

            %acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            %acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (icf & acf) | (tecf & acf);
            exploit = exploit & ~explore;

            % initiation point for action sequences, both exploration and exploitation
            starts = exploit | explore | ps;

            % convert to timestamps
            exploit_ix = find(exploit);
            explore_ix = find(explore);
            starts_ix = find(starts);

            idx = 0;

            idx = idx + 1;
            multi.names{idx} = 'frames';
            multi.onsets{idx} = regs.timestamps';
            multi.durations{idx} = regs.durations';

            multi.orth{idx} = 0; % do not orthogonalise them

            pidx = 0;

            % add exploitation ramps 
            if any(exploit_ix)
                ramp = zeros(size(exploit));
                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ix)
                    eix = exploit_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    ramp(six:eix) = linspace(1, 0, eix - six + 1);
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'exploit';
                multi.pmod(idx).param{pidx} = ramp;
                multi.pmod(idx).poly{pidx} = 1;
            end

            % add exploration ramps 
            if any(explore_ix)
                ramp = zeros(size(explore));
                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ix)
                    eix = explore_ix(i);
                    sixx = find(starts_ix < eix);
                    six = starts_ix(sixx(end));
                    ramp(six:eix) = linspace(1, 0, eix - six + 1);
                end
                pidx = pidx + 1;
                multi.pmod(idx).name{pidx} = 'explore';
                multi.pmod(idx).param{pidx} = ramp;
                multi.pmod(idx).poly{pidx} = 1;
            end


            % optional control regressors
            if ismember(glmodel, [146, 148])
                if any(regs.theory_change_flag_onsets)
                    idx = idx + 1;
                    multi.names{idx} = 'theory_change_flag';
                    multi.onsets{idx} = regs.theory_change_flag_onsets;
                    multi.durations{idx} = zeros(size(multi.onsets{idx}));;
                end
            end
            if ismember(glmodel, [147, 148])
                multi = add_keypresses_to_multi(multi, subj_id, run, conn);
            end

        % exploration versus exploitation, impulses at onsets and offsets
        case {149, 150, 152, 153}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            icf = logical(regs.interaction_change_flag);
            tecf = logical(regs.termination_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps); % cross-reference them; annoying stuff
            assert(all(which));

            %acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            %acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            if ismember(glmodel, [152, 153])
                explore = (tcf & acf);
            else
                explore = (icf & acf) | (tecf & acf);
            end
            exploit = exploit & ~explore; % orthogonalize them

            % initiation point for action sequences, both exploration and exploitation
            starts = exploit | explore | ps;

            % convert to timestamps
            exploit_ts = regs.timestamps(exploit);
            explore_ts = regs.timestamps(explore);
            starts_ts = regs.timestamps(starts);

            idx = 0;

            % add exploitation boxcars
            if any(exploit_ts)
                idx = idx + 1;
                multi.names{idx} = 'exploit_start';
                multi.onsets{idx} = nan(size(exploit_ts));
                multi.durations{idx} = zeros(size(multi.onsets{idx}));

                idx = idx + 1;
                multi.names{idx} = 'exploit_end';
                multi.onsets{idx} = nan(size(exploit_ts));
                multi.durations{idx} = zeros(size(multi.onsets{idx}));

                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ts)
                    en = exploit_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx - 1}(i) = st;
                    multi.onsets{idx}(i) = en;
                end
            end

            % add exploration boxcars
            if any(explore_ts)
                idx = idx + 1;
                multi.names{idx} = 'explore_start';
                multi.onsets{idx} = nan(size(explore_ts));
                multi.durations{idx} = zeros(size(multi.onsets{idx}));

                idx = idx + 1;
                multi.names{idx} = 'explore_end';
                multi.onsets{idx} = nan(size(explore_ts));
                multi.durations{idx} = zeros(size(multi.onsets{idx}));

                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ts)
                    en = explore_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx - 1}(i) = st;
                    multi.onsets{idx}(i) = en;
                end
            end

            % optional control regressors
            if ismember(glmodel, [150, 153])
                % GLM 9: nuisance regressors
                multi = add_games_to_multi(multi, subj_id, run, conn);
                multi = add_keyholds_to_multi(multi, subj_id, run, conn);
                multi = add_visuals_to_multi(multi, subj_id, run, conn);
                multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});
            end


        % exploration versus exploitation, impulses at onsets; # keypresses pmods; subgoal completion
        case {151}

            [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn, true);

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            tcf = logical(regs.theory_change_flag);
            icf = logical(regs.interaction_change_flag);
            tecf = logical(regs.termination_change_flag);
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps); % cross-reference them; annoying stuff
            assert(all(which));

            %acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            %acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (icf & acf) | (tecf & acf);
            exploit = exploit & ~explore; % orthogonalize them

            % initiation point for action sequences, both exploration and exploitation
            starts = exploit | explore | ps;
            ends = exploit | explore;

            % convert to timestamps
            exploit_ts = regs.timestamps(exploit);
            explore_ts = regs.timestamps(explore);
            starts_ts = regs.timestamps(starts);
            ends_ts = regs.timestamps(ends);

            idx = 0;

            % add exploitation impulses
            if any(exploit_ts)
                idx = idx + 1;
                multi.names{idx} = 'exploit_start';
                multi.onsets{idx} = nan(size(exploit_ts));
                multi.durations{idx} = zeros(size(multi.onsets{idx}));
                num_keypresses = zeros(size(multi.onsets{idx}));

                multi.orth{idx} = 0; % do not orthogonalise them

                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ts)
                    en = exploit_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx}(i) = st;

                    % Parametrically modulated by number of key presses until the subgoal
                    for k = 1:length(keypresses)
                        num_keypresses(i) = num_keypresses(i) + sum((keypresses{k} >= st) & (keypresses{k} <= en));
                    end
                end

                if std(num_keypresses) > 1e-6
                    % Make sure they're not the same 
                    multi.pmod(idx).name{1} = 'keypresses';
                    multi.pmod(idx).param{1} = num_keypresses;
                    multi.pmod(idx).poly{1} = 1;
                end
            end

            % add exploration impulses
            if any(explore_ts)
                idx = idx + 1;
                multi.names{idx} = 'explore_start';
                multi.onsets{idx} = nan(size(explore_ts));
                multi.durations{idx} = zeros(size(multi.onsets{idx}));
                num_keypresses = zeros(size(multi.onsets{idx}));

                multi.orth{idx} = 0; % do not orthogonalise them

                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ts)
                    en = explore_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx}(i) = st;

                    % Parametrically modulated by number of key presses until the subgoal
                    for k = 1:length(keypresses)
                        num_keypresses(i) = num_keypresses(i) + sum((keypresses{k} >= st) & (keypresses{k} <= en));
                    end
                end

                if std(num_keypresses) > 1e-6
                    % Make sure they're not the same 
                    multi.pmod(idx).name{1} = 'keypresses';
                    multi.pmod(idx).param{1} = num_keypresses;
                    multi.pmod(idx).poly{1} = 1;
                end
            end

            % add subgoal completion impulses
            if any(ends_ts)
                idx = idx + 1;
                multi.names{idx} = 'subgoal';
                multi.onsets{idx} = ends_ts;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));
            end

            % add theories updates control regressor
            if any(regs.theory_change_flag_onsets)
                idx = idx + 1;
                multi.names{idx} = 'theory_change_flag';
                multi.onsets{idx} = regs.theory_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


        % hypothesized_terminations_flag
        %
        case 154

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            ht = logical(regs.hypothesized_terminations_flag);
            hypothesized_terminations_flag_onsets = regs.timestamps(ht);

            multi.names{1} = 'hypothesized_terminations_flag';
            multi.onsets{1} = hypothesized_terminations_flag_onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));

        % falsified_terminations_flag
        %
        case 155

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            ft = logical(regs.falsified_terminations_flag);
            falsified_terminations_flag_onsets = regs.timestamps(ft);

            multi.names{1} = 'falsified_terminations_flag';
            multi.onsets{1} = falsified_terminations_flag_onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));

        % new_termination_change_flag = hypothesized_terminations_flag | falsified_terminations_flag
        %
        case 156

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            ht = logical(regs.hypothesized_terminations_flag);
            ft = logical(regs.falsified_terminations_flag);
            new_termination_change_flag = ht | ft;
            new_termination_change_flag_onsets = regs.timestamps(new_termination_change_flag);

            multi.names{1} = 'new_termination_change_flag';
            multi.onsets{1} = new_termination_change_flag_onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));

        % GLM 102, with new_termination_change_flag = hypothesized_terminations_flag | falsified_terminations_flag
        % and theory_change_flag = sprite_change_flag | interaction_change_flag | new_termination_change_flag
        %
        case 157

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            sc = logical(regs.sprite_change_flag);
            ic = logical(regs.interaction_change_flag);
            ht = logical(regs.hypothesized_terminations_flag);
            ft = logical(regs.falsified_terminations_flag);
            new_theory_change_flag = sc | ic | ht | ft;
            new_theory_change_flag_onsets = regs.timestamps(new_theory_change_flag);

            multi.names{1} = 'new_theory_change_flag';
            multi.onsets{1} = new_theory_change_flag_onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 102, with 
        % and theory_change_flag = sprite_change_flag | interaction_change_flag | falsified_terminations  (i.e. what Pedro thinks it should be)
        %
        case 158

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            sc = logical(regs.sprite_change_flag);
            ic = logical(regs.interaction_change_flag);
            ft = logical(regs.falsified_terminations_flag);
            new_theory_change_flag = sc | ic | ft;
            new_theory_change_flag_onsets = regs.timestamps(new_theory_change_flag);

            multi.names{1} = 'new_theory_change_flag';
            multi.onsets{1} = new_theory_change_flag_onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});


        % GLM 103, with empa_regressors
        %
        case 159

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');

            if length(regs.sprite_change_flag_onsets) > 0
                multi.names{1} = 'sprite_change_flag';
                multi.onsets{1} = regs.sprite_change_flag_onsets;
                multi.durations{1} = zeros(size(multi.onsets{1}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 104, with empa_regressors
        %
        case 160

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');

            if length(regs.interaction_change_flag_onsets) > 0
                multi.names{1} = 'interaction_change_flag';
                multi.onsets{1} = regs.interaction_change_flag_onsets;
                multi.durations{1} = zeros(size(multi.onsets{1}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});


        % GLM 105, with new_termination_change_flag = hypothesized_terminations_flag | falsified_terminations_flag
        %
        case 161

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            ht = logical(regs.hypothesized_terminations_flag);
            ft = logical(regs.falsified_terminations_flag);
            new_termination_change_flag = ht | ft;
            new_termination_change_flag_onsets = regs.timestamps(new_termination_change_flag);

            if length(new_termination_change_flag_onsets) > 0
                multi.names{1} = 'new_termination_change_flag';
                multi.onsets{1} = new_termination_change_flag_onsets;
                multi.durations{1} = zeros(size(multi.onsets{1}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 105, falsified_terminations_flag
        %
        case 162

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            ft = logical(regs.falsified_terminations_flag);
            falsified_terminations_flag_onsets = regs.timestamps(ft);

            if length(falsified_terminations_flag_onsets) > 0
                multi.names{1} = 'falsified_terminations_flag';
                multi.onsets{1} = falsified_terminations_flag_onsets;
                multi.durations{1} = zeros(size(multi.onsets{1}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 106, but with new_termination_change_flag = hypothesized_terminations_flag | falsified_terminations_flag
        case 163

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            ht = logical(regs.hypothesized_terminations_flag);
            ft = logical(regs.falsified_terminations_flag);
            new_termination_change_flag = ht | ft;
            new_termination_change_flag_onsets = regs.timestamps(new_termination_change_flag);
            if length(new_termination_change_flag_onsets) > 0 
                idx = idx + 1;
                multi.names{idx} = 'new_termination_change_flag';
                multi.onsets{idx} = new_termination_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % GLM 106, but with falsified_terminations_flag
        case 164

            idx = 0;

            regs = get_regressors(subj_id, run, conn, true, 'empa_regressors');

            if length(regs.sprite_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'sprite_change_flag';
                multi.onsets{idx} = regs.sprite_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            if length(regs.interaction_change_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'interaction_change_flag';
                multi.onsets{idx} = regs.interaction_change_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end


            ht = logical(regs.hypothesized_terminations_flag);
            ft = logical(regs.falsified_terminations_flag);
            falsified_terminations_flag_onsets = regs.timestamps(ft);
            if length(falsified_terminations_flag_onsets) > 0
                idx = idx + 1;
                multi.names{idx} = 'falsified_terminations_flag';
                multi.onsets{idx} = falsified_terminations_flag_onsets;
                multi.durations{idx} = zeros(size(multi.onsets{idx}));;
            end

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % planning @ avatar collision flag: win_plan_length
        %
        case {165, 166}
            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            win_plan_length = regs.win_plan_length;
            win_plan_length(isnan(win_plan_length)) = 0;

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);

            acf = acf | (win_plan_length > 0); % include beginning-of-episode planning

            multi.names{1} = 'avatar_collision_flag';
            multi.onsets{1} = regs.timestamps(acf);
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            multi.orth{1} = 0; % do not orthogonalise them

            switch glmodel
                case 165
                    multi.pmod(1).name{1} = 'win_plan_length';
                    multi.pmod(1).param{1} = win_plan_length(acf);
                case 166
                    multi.pmod(1).name{1} = 'log_win_plan_length';
                    multi.pmod(1).param{1} = log(win_plan_length(acf));
                    % cleanup nans/infs
                    multi.pmod(1).param{1}(isnan(multi.pmod(1).param{1})) = 0;
                    multi.pmod(1).param{1}(isinf(multi.pmod(1).param{1})) = 0;
                otherwise
                    assert(false);
            end
            multi.pmod(1).poly{1} = 1;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % planning @ avatar collision flag: avg_plan_length
        %
        case {167, 168}
            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            avg_plan_length = regs.avg_plan_length;
            avg_plan_length(isnan(avg_plan_length)) = 0;

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);

            acf = acf | (avg_plan_length > 0); % include beginning-of-episode planning

            multi.names{1} = 'avatar_collision_flag';
            multi.onsets{1} = regs.timestamps(acf);
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            multi.orth{1} = 0; % do not orthogonalise them

            switch glmodel
                case 167
                    multi.pmod(1).name{1} = 'avg_plan_length';
                    multi.pmod(1).param{1} = avg_plan_length(acf);
                case 168
                    multi.pmod(1).name{1} = 'log_avg_plan_length';
                    multi.pmod(1).param{1} = log(avg_plan_length(acf));
                    % cleanup nans/infs
                    multi.pmod(1).param{1}(isnan(multi.pmod(1).param{1})) = 0;
                    multi.pmod(1).param{1}(isinf(multi.pmod(1).param{1})) = 0;
                otherwise
                    assert(false);
            end
            multi.pmod(1).poly{1} = 1;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % planning @ avatar collision flag: num_plans
        %
        case {169, 170}
            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            num_plans = regs.num_plans;
            num_plans(isnan(num_plans)) = 0;

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);

            acf = acf | (num_plans > 0); % include beginning-of-episode planning

            multi.names{1} = 'avatar_collision_flag';
            multi.onsets{1} = regs.timestamps(acf);
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            multi.orth{1} = 0; % do not orthogonalise them

            switch glmodel
                case 169
                    multi.pmod(1).name{1} = 'num_plans';
                    multi.pmod(1).param{1} = num_plans(acf);
                case 170
                    multi.pmod(1).name{1} = 'log_num_plans';
                    multi.pmod(1).param{1} = log(num_plans(acf));
                    % cleanup nans/infs
                    multi.pmod(1).param{1}(isnan(multi.pmod(1).param{1})) = 0;
                    multi.pmod(1).param{1}(isinf(multi.pmod(1).param{1})) = 0;
                otherwise
                    assert(false);
            end
            multi.pmod(1).poly{1} = 1;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % planning @ avatar collision flag: total_plan_length = avg_plan_length * num_plans
        %
        case {171, 172}
            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            num_plans = regs.num_plans;
            avg_plan_length = regs.avg_plan_length;
            num_plans(isnan(num_plans)) = 0;
            avg_plan_length(isnan(avg_plan_length)) = 0;

            total_plan_length = num_plans .* avg_plan_length;

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);

            acf = acf | (total_plan_length > 0); % include beginning-of-episode planning

            multi.names{1} = 'avatar_collision_flag';
            multi.onsets{1} = regs.timestamps(acf);
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            multi.orth{1} = 0; % do not orthogonalise them

            switch glmodel
                case 171
                    multi.pmod(1).name{1} = 'total_plan_length';
                    multi.pmod(1).param{1} = total_plan_length(acf);
                case 172
                    multi.pmod(1).name{1} = 'log_total_plan_length';
                    multi.pmod(1).param{1} = log(total_plan_length(acf));
                    % cleanup nans/infs
                    multi.pmod(1).param{1}(isnan(multi.pmod(1).param{1})) = 0;
                    multi.pmod(1).param{1}(isinf(multi.pmod(1).param{1})) = 0;
                otherwise
                    assert(false);
            end
            multi.pmod(1).poly{1} = 1;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % planning @ avatar collision flag: win_plan_ac
        %
        case {173, 174}
            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            win_plan_ac = regs.win_plan_ac;
            win_plan_ac(isnan(win_plan_ac)) = 0;

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);

            acf = acf | (win_plan_ac > 0); % include beginning-of-episode planning

            multi.names{1} = 'avatar_collision_flag';
            multi.onsets{1} = regs.timestamps(acf);
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            multi.orth{1} = 0; % do not orthogonalise them

            switch glmodel
                case 173
                    multi.pmod(1).name{1} = 'win_plan_ac';
                    multi.pmod(1).param{1} = win_plan_ac(acf);
                case 174
                    multi.pmod(1).name{1} = 'log_win_plan_ac';
                    multi.pmod(1).param{1} = log(win_plan_ac(acf));
                    % cleanup nans/infs
                    multi.pmod(1).param{1}(isnan(multi.pmod(1).param{1})) = 0;
                    multi.pmod(1).param{1}(isinf(multi.pmod(1).param{1})) = 0;
                otherwise
                    assert(false);
            end
            multi.pmod(1).poly{1} = 1;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % planning @ avatar collision flag: avg_plan_ac
        %
        case {175, 176}
            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            avg_plan_ac = regs.avg_plan_ac;
            avg_plan_ac(isnan(avg_plan_ac)) = 0;

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);

            acf = acf | (avg_plan_ac > 0); % include beginning-of-episode planning

            multi.names{1} = 'avatar_collision_flag';
            multi.onsets{1} = regs.timestamps(acf);
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            multi.orth{1} = 0; % do not orthogonalise them

            switch glmodel
                case 175
                    multi.pmod(1).name{1} = 'avg_plan_ac';
                    multi.pmod(1).param{1} = avg_plan_ac(acf);
                case 176
                    multi.pmod(1).name{1} = 'log_avg_plan_ac';
                    multi.pmod(1).param{1} = log(avg_plan_ac(acf));
                    % cleanup nans/infs
                    multi.pmod(1).param{1}(isnan(multi.pmod(1).param{1})) = 0;
                    multi.pmod(1).param{1}(isinf(multi.pmod(1).param{1})) = 0;
                otherwise
                    assert(false);
            end
            multi.pmod(1).poly{1} = 1;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % planning @ avatar collision flag: win_plan_eff
        %
        case {177, 178}
            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            win_plan_eff = regs.win_plan_eff;
            win_plan_eff(isnan(win_plan_eff)) = 0;

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);

            acf = acf | (win_plan_eff > 0); % include beginning-of-episode planning

            multi.names{1} = 'avatar_collision_flag';
            multi.onsets{1} = regs.timestamps(acf);
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            multi.orth{1} = 0; % do not orthogonalise them

            switch glmodel
                case 177
                    multi.pmod(1).name{1} = 'win_plan_eff';
                    multi.pmod(1).param{1} = win_plan_eff(acf);
                case 178
                    multi.pmod(1).name{1} = 'log_win_plan_eff';
                    multi.pmod(1).param{1} = log(win_plan_eff(acf));
                    % cleanup nans/infs
                    multi.pmod(1).param{1}(isnan(multi.pmod(1).param{1})) = 0;
                    multi.pmod(1).param{1}(isinf(multi.pmod(1).param{1})) = 0;
                otherwise
                    assert(false);
            end
            multi.pmod(1).poly{1} = 1;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % planning @ avatar collision flag: avg_plan_eff
        %
        case {179, 180}
            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            avg_plan_eff = regs.avg_plan_eff;
            avg_plan_eff(isnan(avg_plan_eff)) = 0;

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);

            acf = acf | (avg_plan_eff > 0); % include beginning-of-episode planning

            multi.names{1} = 'avatar_collision_flag';
            multi.onsets{1} = regs.timestamps(acf);
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            multi.orth{1} = 0; % do not orthogonalise them

            switch glmodel
                case 179
                    multi.pmod(1).name{1} = 'avg_plan_eff';
                    multi.pmod(1).param{1} = avg_plan_eff(acf);
                case 180
                    multi.pmod(1).name{1} = 'log_avg_plan_eff';
                    multi.pmod(1).param{1} = log(avg_plan_eff(acf));
                    % cleanup nans/infs
                    multi.pmod(1).param{1}(isnan(multi.pmod(1).param{1})) = 0;
                    multi.pmod(1).param{1}(isinf(multi.pmod(1).param{1})) = 0;
                otherwise
                    assert(false);
            end
            multi.pmod(1).poly{1} = 1;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % exploration versus exploitation with dTnov, boxcars
        % copy of 116 et al
        case {181, 182, 183, 184}

            [regs, ~, ~] = get_regressors(subj_id, run, conn, true, 'empa_regressors');
            dTnov_len = regs.dTnov_len;
            sgf1 = logical(regs.subgoal_flag1);
            sgf2 = logical(regs.subgoal_flag2);

            [~, visuals, ~] = get_visuals(subj_id, run, conn, true, 'empa_plays_post3');
            acf = logical(visuals.avatar_collision_flag);
            ks = logical(visuals.killed_sprites);
            ps = logical(visuals.play_start);
            pe = logical(visuals.play_end);
            win = logical(visuals.win);

            which = ismember(visuals.timestamps, regs.state_timestamps); % cross-reference them; annoying stuff
            if ~all(which)
                % TODO this seems outdated, never hits; remove eventually
                keyboard
            end
            acf = acf(which);
            ks = ks(which);
            ps = ps(which);
            pe = pe(which);
            win = win(which);

            %acf_broad = acf | [0; acf(1:end-1)];  % account for off-by-one
            %acf = acf_broad;

            % exploration = subgoals that bring us closer to satisfying termination conditions
            % exploration = subgoals that lead to theory updating
            exploit = (sgf1 & acf) | (sgf2 & acf & ks) | (win & acf);
            explore = (dTnov_len < 0);
            exploit = exploit & ~explore; % orthogonalize them

            % initiation point for action sequences, both exploration and exploitation
            switch glmodel
                case {181, 183}
                    starts = exploit | explore | ps;
                case {182, 184}
                    starts = exploit | explore | ps | acf;
                otherwise
                    assert(false);
            end

            % convert to timestamps
            exploit_ts = regs.timestamps(exploit);
            explore_ts = regs.timestamps(explore);
            starts_ts = regs.timestamps(starts);

            idx = 0;

            % add exploitation boxcars
            if any(exploit_ts)
                idx = idx + 1;
                multi.names{idx} = 'exploit';
                multi.onsets{idx} = nan(size(exploit_ts));
                multi.durations{idx} = nan(size(multi.onsets{idx}));

                % match exploitation subgoal with closest proceeding action initiation
                for i = 1:length(exploit_ts)
                    en = exploit_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx}(i) = st;
                    multi.durations{idx}(i) = en - st;
                end
            end

            % add exploration boxcars
            if any(explore_ts)
                idx = idx + 1;
                multi.names{idx} = 'explore';
                multi.onsets{idx} = nan(size(explore_ts));
                multi.durations{idx} = nan(size(multi.onsets{idx}));

                % match exploration subgoal with closest proceeding action initiation
                for i = 1:length(explore_ts)
                    en = explore_ts(i);
                    ix = find(starts_ts < en);
                    st = starts_ts(ix(end));
                    multi.onsets{idx}(i) = st;
                    multi.durations{idx}(i) = en - st;
                end
            end

            % optional control regressors
            switch glmodel
                case {181, 182}
                    % GLM 9: nuisance regressors
                    multi = add_games_to_multi(multi, subj_id, run, conn);
                    multi = add_keyholds_to_multi(multi, subj_id, run, conn);
                    multi = add_visuals_to_multi(multi, subj_id, run, conn);
                    multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

                case {183, 184}
                    if any(regs.theory_change_flag_onsets)
                        idx = idx + 1;
                        multi.names{idx} = 'theory_change_flag';
                        multi.onsets{idx} = regs.theory_change_flag_onsets;
                        multi.durations{idx} = zeros(size(multi.onsets{idx}));;
                    end

                otherwise
                    assert(false);
            end

        % like GLM 102, but with lags
        %
        case {185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195}

            % optionally offset from the theory_change_flag time
            lag = (glmodel - 185) * 2;
            fprintf('glm %d -> lag %d', glmodel, lag);

            % from GLM 3: theory_change_flag
            %
            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.theory_change_flag_onsets;

            multi.names{1} = 'theory_change_flag';
            multi.onsets{1} = onsets + lag;
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        % like GLM 102, but with separate regressors for each data partition
        %
        case 196

            % from GLM 3: theory_change_flag
            %
            regs = get_regressors(subj_id, run, conn, true);
            onsets = regs.theory_change_flag_onsets;

            partition_id = partition_id_from_run_id(run_id);
            multi.names{1} = sprintf('theory_change_flag_part%d', partition_id);
            multi.onsets{1} = onsets;
            multi.durations{1} = zeros(size(multi.onsets{1}));;

            % GLM 9: nuisance regressors
            multi = add_games_to_multi(multi, subj_id, run, conn);
            multi = add_keyholds_to_multi(multi, subj_id, run, conn);
            multi = add_visuals_to_multi(multi, subj_id, run, conn);
            multi = add_onoff_to_multi(multi, subj_id, run, conn, {'play_start', 'play_end'});

        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

    if save_output
        save(filename, 'multi', '-v7.3'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
    end

    if exist('conn', 'var')
        close(conn);
    end

end

% helper functions that augment the multi structure 

function multi = add_games_to_multi(multi, subj_id, run, conn)
    % from GLM 1: game instance boxcar regressor
    %
    if isfield(multi, 'names')
        assert(length(multi.names) == length(multi.onsets));
        assert(length(multi.names) == length(multi.durations));
        idx = length(multi.names);
    else
        idx = 0;
    end

    [game_names, onsets, durs] = get_games(subj_id, run, conn);
    for i = 1:numel(game_names)
        % instance boxcar regressor
        idx = idx + 1;
        multi.names{idx} = game_names{i};
        multi.onsets{idx} = onsets{i};
        multi.durations{idx} = durs{i};
    end
end

function multi = add_keypresses_to_multi(multi, subj_id, run, conn)
    if isfield(multi, 'names')
        assert(length(multi.names) == length(multi.onsets));
        assert(length(multi.names) == length(multi.durations));
        idx = length(multi.names);
    else
        idx = 0;
    end

    [keyNames, keyholds, keyholds_post, keypresses] = get_keypresses(subj_id, run, conn, true);

    for k = 1:numel(keyNames)
        if size(keypresses{k}, 1) > 0
            idx = idx + 1;
            multi.names{idx} = keyNames{k}; % TODO rm prefix... but already ran
            multi.onsets{idx} = keypresses{k}(:,1)';
            multi.durations{idx} = zeros(size(multi.onsets{idx}));
        end
    end
end

function multi = add_keyholds_to_multi(multi, subj_id, run, conn)
    % from GLM 5: keyholds nuisance regressors
    %
    if isfield(multi, 'names')
        assert(length(multi.names) == length(multi.onsets));
        assert(length(multi.names) == length(multi.durations));
        idx = length(multi.names);
    else
        idx = 0;
    end

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
end

function multi = add_visuals_to_multi(multi, subj_id, run, conn, fields)
    % GLM 7: frame nuisance regressors
    %
    if isfield(multi, 'names')
        assert(length(multi.names) == length(multi.onsets));
        assert(length(multi.names) == length(multi.durations));
        idx = length(multi.names);
    else
        idx = 0;
    end

    % Notice that we use legacy_fields for backwards compatibility / reproducibility with all the analyses
    [legacy_fields, visuals] = get_visuals(subj_id, run, conn, true);

    if ~exist('fields', 'var')
        fields = legacy_fields;
    end

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
            if ((subj_id == 18 && run.run_id == 2) || (subj_id == 30 && run.run_id == 4)) && strcmp(fields{i}, 'effectsByCol')
                % collinearity, special case
                continue
            end
            pix = pix + 1;
            multi.pmod(idx).name{pix} = fields{i};
            multi.pmod(idx).param{pix} = visuals.(fields{i});
            multi.pmod(idx).poly{pix} = 1;
        end
    end
end

function multi = add_onoff_to_multi(multi, subj_id, run, conn, fields)
    % GLM 8: on/off nuisance regressors
    %
    if isfield(multi, 'names')
        assert(length(multi.names) == length(multi.onsets));
        assert(length(multi.names) == length(multi.durations));
        idx = length(multi.names);
    else
        idx = 0;
    end

    [onoff] = get_onoff(subj_id, run, conn, true);
   
    if ~exist('fields', 'var')
        fields = fieldnames(onoff);
    end

    for i = 1:numel(fields)
        idx = idx + 1;
        multi.names{idx} = fields{i};
        multi.onsets{idx} = onoff.(fields{i});
        multi.durations{idx} = zeros(size(multi.onsets{idx}));
    end
end

function multi = add_frames_with_pmods(multi, subj_id, run, conn, fields)
    % add regressors as parametric modulators to frames
    % works for regressors with a single variable (e.g. not R_GGs)

    [regs, ~, ~] = get_regressors(subj_id, run, conn, true);

    if isfield(multi, 'names')
        idx = length(multi.names) + 1;
    else
        idx = 1;
    end
    multi.names{idx} = 'frames';
    multi.onsets{idx} = regs.timestamps';
    multi.durations{idx} = regs.durations';

    multi.orth{idx} = 0; % do not orthogonalise them

    for i = 1:numel(fields)
        field = fields{i};

        assert(size(regs.(field), 2) == 1);
        if strcmp(field, 'spriteKL') && max(regs.(field)) - min(regs.(field)) < 1e-3
            % all values are basically the same (0)
            % if we check for strict equality, will get invalid contrasts
            continue;
        end

        multi.pmod(idx).name{i} = field;
        multi.pmod(idx).param{i} = regs.(field);
        multi.pmod(idx).poly{i} = 1;

        multi.pmod(idx).param{i}(isnan(multi.pmod(idx).param{i})) = 0; % TODO happens for lik during first few timesteps; ideally, remove altogether
    end
end
