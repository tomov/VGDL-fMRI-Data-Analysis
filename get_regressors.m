function [regs, X, fields] = get_regressors(subj_id, run, conn, do_cache, regressors_collection)

    %{
    subj_id = 1; 
    run_id = 2; 
    conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54');
    query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id); 
    run = find(conn, 'runs', 'query', query);
    regressors_collection = 'regressors';
    do_cache = false;
    %}

    if ~exist('do_cache', 'var')
        do_cache = false;
    end

    if ~exist('regressors_collection', 'var')
        regressors_collection = 'regressors';
        filename = fullfile(get_mat_dir(false), sprintf('get_regressors_subj%d_run%d.mat', subj_id, run.run_id));
    else
        filename = fullfile(get_mat_dir(false), sprintf('get_regressors_subj%d_run%d_c=%s.mat', subj_id, run.run_id, regressors_collection));
    end
    filename

    % optionally cache
    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    % helper function to get regressors for each frame in vgdl_create_multi
    % copied & modified from GLM 3
    % note run is a struct
    %

    assert(isequal(run.subj_id, num2str(subj_id)));

    plays_post_collection = 'plays_post'; % by default, unless specified otherwise
    switch regressors_collection
        case {'empa_regressors'}
            % 818de8d483c140e54eee1c57fae0c0fc1c457725
            plays_post_collection = 'empa_plays_post2'; % separate collection

            % TODO termination_change_flag & interaction_change_flag differ between regressors and plays_post; latter seem more correct
            binreg_fields = {'theory_change_flag', 'sprite_change_flag', 'newEffects_flag', 'newTimeStep_flag', 'replan_flag'}; % binary db.regressors => onsets only, durations irrelevant; have the option of having them as onsets only
            reg_fields = [binreg_fields, {'spriteKL', 'sum_lik', 'n_ts', 'num_effects', 'R_GG', 'R_GGs', 'R_SG', 'R_SGs'}]; % db.regressors 
            binpost_fields = {'interaction_change_flag', 'termination_change_flag'}; % binary db.plays_post, copied & fixed from db.regressors
            post_fields = [binpost_fields {'likelihood', 'surprise', 'sum_lik_play', 'S_len','I_len','T_len','Igen_len','Tnov_len','Ip_len','dS_len','dI_len','dT_len','dIgen_len','dTnov_len','dIp_len', 'subgoal_flag1', 'subgoal_flag2'}]; % db.plays_post 

            % empa_plays_post2
            post_fields = [post_fields {'hypothesized_new_terminations_len', 'hypothesized_new_terminations_flag', 'falsified_existing_terminations_len', 'falsified_existing_terminations_flag'}];

        case {'regressors', 'regressors_DELETEME'}
            % current one -- regressors_and_playspost_2020_05_31_finalTS_block

            % TODO termination_change_flag & interaction_change_flag differ between regressors and plays_post; latter seem more correct
            binreg_fields = {'theory_change_flag', 'sprite_change_flag', 'newEffects_flag', 'newTimeStep_flag', 'replan_flag'}; % binary db.regressors => onsets only, durations irrelevant; have the option of having them as onsets only
            reg_fields = [binreg_fields, {'spriteKL', 'sum_lik', 'n_ts', 'num_effects', 'R_GG', 'R_GGs', 'R_SG', 'R_SGs'}]; % db.regressors 
            binpost_fields = {'interaction_change_flag', 'termination_change_flag'}; % binary db.plays_post, copied & fixed from db.regressors
            post_fields = [binpost_fields {'likelihood', 'surprise', 'sum_lik_play', 'S_len','I_len','T_len','Igen_len','Tnov_len','Ip_len','dS_len','dI_len','dT_len','dIgen_len','dTnov_len','dIp_len'}]; % db.plays_post 

        case 'regressors_1'

            % oldest one -- regressors_2020_04_01_fullTSList_nonrefactor

            binreg_fields = {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};
            reg_fields = [binreg_fields, {}]  
            binpost_fields = {};
            post_fields = [binpost_fields {}]; % db.plays_post 

        case 'regressors_2'

            %  regressors_2020_05_19_finalTimeStep_reset_nonrefactored

            binreg_fields = {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};
            reg_fields = [binreg_fields, {}]  
            binpost_fields = {};
            post_fields = [binpost_fields {}]; % db.plays_post 


        case 'regressors_3'

            %  regressors_2020_05_24_finalTimeStep_reset_refactored

            binreg_fields = {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};
            reg_fields = [binreg_fields, {}]  
            binpost_fields = {};
            post_fields = [binpost_fields {}]; % db.plays_post 

        % TODO create regressors_4 from regressors_and_playspost_2020_05_31_finalTS_reset2_refactored

        otherwise

            assert(false);
    end

    
    %{
    % for loading old regressors
    binreg_fields = {'theory_change_flag'};
    reg_fields = {};
    binpost_fields = {};
    post_fields = {};
    %}

    % TODO dedupe w/ vgdl_create_rsa and others
    game_names_ordered = get_game_names_ordered(subj_id);
    game_name_to_id = containers.Map(game_names_ordered, 1:6);

    regs = struct;
    for i = 1:numel(binreg_fields)
        regs.([binreg_fields{i}, '_onsets']) = [];
    end
    for i = 1:numel(reg_fields)
        regs.(reg_fields{i}) = [];
    end
    for i = 1:numel(binpost_fields)
        regs.([binpost_fields{i}, '_onsets']) = [];
    end
    for i = 1:numel(post_fields)
        regs.(post_fields{i}) = [];
    end
    regs.durations = [];
    regs.game_ids = [];
    regs.block_ids = [];
    regs.instance_ids = [];
    regs.play_ids = [];
    regs.keystate_timestamps = []; % these all come from EMPA, so the timestamps are keystate['ts'] timestamps, whereas those from the states have state['ts'] timestamps; TODO unify maybe
    regs.state_timestamps = []; % these all come directly from states, so they are state['ts'] timestamps;

    regs.game_names_ordered = game_names_ordered;
    regs.game_name_to_id = game_name_to_id;

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

                if ismember(regressors_collection, {'regressors', 'empa_regressors'})
                    regressors = find(conn, regressors_collection, 'query', q);
                    assert(length(regressors) == 1);
                else
                    % we didn't have our shit together back in the day
                    regressors = find(conn, regressors_collection, 'query', q, 'sort', '{"dt": -1.0}'); % momchil: assume latest one is the correct one 
                    if isempty(regressors)
                        disp('EMPTY w t f........');
                        continue
                    end
                end
                reg = regressors(1);
   
                % some plays are too short (e.g. 1-2 frames)
                if length(reg.regressors.theory_change_flag) == 0
                    continue
                end

                assert(immse(cellfun(@(x) x{3}, reg.regressors.theory_change_flag), cellfun(@(x) x{3}, reg.regressors.interaction_change_flag)) < 1e-10); % assert identical timestamps
                assert(immse(cellfun(@(x) x{3}, reg.regressors.theory_change_flag), cellfun(@(x) x{3}, reg.regressors.sprite_change_flag)) < 1e-10); % assert identical timestamps
                assert(immse(cellfun(@(x) x{3}, reg.regressors.theory_change_flag), cellfun(@(x) x{3}, reg.regressors.termination_change_flag)) < 1e-10); % assert identical timestamps

                for i = 1:numel(binreg_fields)
                    if iscell(reg.regressors.(binreg_fields{i}))
                        r = cellfun(@(x) x{1}, reg.regressors.(binreg_fields{i}));
                        t = cellfun(@(x) x{3}, reg.regressors.(binreg_fields{i}));
                    else
                        r = reg.regressors.(binreg_fields{i})(:,1);
                        t = reg.regressors.(binreg_fields{i})(:,3);
                    end
                    regs.([binreg_fields{i}, '_onsets']) = [regs.([binreg_fields{i}, '_onsets']); t(find(r > 0)) - run.scan_start_ts];
                end

                for i = 1:numel(reg_fields)
                    if iscell(reg.regressors.(reg_fields{i}))
                        if numel(reg.regressors.(reg_fields{i}){1}{1}) == 1
                            r = cellfun(@(x) x{1}, reg.regressors.(reg_fields{i}));
                        else
                            r = cellfun(@(x) x{1}', reg.regressors.(reg_fields{i}), 'UniformOutput', false);
                            r = cell2mat(r);
                        end
                    else
                        r = reg.regressors.(reg_fields{i})(:,1);
                    end
                    regs.(reg_fields{i}) = [regs.(reg_fields{i}); r];
                end

                durations = t(2:end) - t(1:end-1);
                durations = [durations; mean(durations)]; % guesstimate duration of last frame
                regs.durations = [regs.durations; durations];
                regs.game_ids = [regs.game_ids; ones(size(durations)) * game_name_to_id(game_name)]; 
                regs.block_ids = [regs.block_ids; ones(size(durations)) * block.block_id]; 
                regs.instance_ids = [regs.instance_ids; ones(size(durations)) * instance.instance_id];
                regs.play_ids = [regs.play_ids; ones(size(durations)) * play.play_id];
                regs.keystate_timestamps = [regs.keystate_timestamps; t - run.scan_start_ts]; % assuming all have the same ts; notice those are keystate timestamps (slightly off from state timestamps... sorry; see core.py)

                ttt = t;

                plays_post = find(conn, plays_post_collection, 'query', q);
                if length(plays_post) ~= 1
                    disp('bbb')
                    keyboard
                end
                assert(length(plays_post) == 1);
                play_post = plays_post(1);

                % assert identical stuff computed during replay and after; note will not be true for older replays 
                %{
                assert(immse(int32(cellfun(@(x) x{1}, reg.regressors.interaction_change_flag)), int32(cellfun(@(x) x{1}, play_post.interaction_change_flag))) < 1e-10); % not equal when interaction set gets reshuffled, probably due to pickling? or mismatch between self.hypotheses[0] and prev_theory in fmri_playsPostproc.py (it changes somewhere?)
                assert(immse(int32(cellfun(@(x) x{1}, reg.regressors.termination_change_flag)), int32(cellfun(@(x) x{1}, play_post.termination_change_flag))) < 1e-10); % same; actually seems like the playsPostproc versions are more legit # TODO look into it!!!
                %}
                %assert(immse(int32(cellfun(@(x) x{1}, reg.regressors.newTimeStep_flag)), int32(play_post.newTimeStep_flag)) < 1e-10);
                %assert(sqrt(nanmean((reg.regressors.likelihood(:,1) - double(play_post.likelihood)).^2)) < 1e-10); % initial NaNs
                %assert(sqrt(nanmean((reg.regressors.surprise(:,1) - double(play_post.surprise)).^2)) < 1e-10); % initial NaNs

                for i = 1:numel(binpost_fields)
                    if iscell(play_post.(binpost_fields{i}))
                        r = cellfun(@(x) x{1}, play_post.(binpost_fields{i}));
                        t = cellfun(@(x) x{3}, play_post.(binpost_fields{i}));
                    else
                        r = play_post.(binpost_fields{i})(:,1);
                        t = play_post.(binpost_fields{i})(:,3);
                    end
                    regs.([binpost_fields{i}, '_onsets']) = [regs.([binpost_fields{i}, '_onsets']); t(find(r > 0)) - run.scan_start_ts];
                end

                for i = 1:numel(post_fields)
                    if iscell(play_post.(post_fields{i}))
                        if numel(play_post.(post_fields{i}){1}{1}) == 1
                            r = cellfun(@(x) x{1}, play_post.(post_fields{i}));
                        else
                            r = cellfun(@(x) x{1}', play_post.(post_fields{i}), 'UniformOutput', false);
                            r = cell2mat(r);
                        end
                    else
                        r = play_post.(post_fields{i})(:,1);
                    end
                    if ismember(post_fields{i}, {'subgoal_flag1', 'subgoal_flag2'})
                        r = r(2:end-1); % these have extra state, just like in get_visuals()
                    end
                    regs.(post_fields{i}) = [regs.(post_fields{i}); r];
                end

                % note we skip first and last state timestamp b/c EMPA skips them too when fMRI logging TODO 
                %if length(play_post.timestamps(2:end-1)) ~= length(t)

                regs.state_timestamps = [regs.state_timestamps; play_post.timestamps(2:end-1) - run.scan_start_ts]; % assuming all have the same ts; notice those are keystate timestamps (slightly off from state timestamps... sorry; see core.py)
                fprintf('  after  %d %d -- %d %d\n', length(regs.state_timestamps), length(regs.keystate_timestamps), length(play_post.timestamps(2:end-1)), length(ttt - run.scan_start_ts)); 

                if length(play_post.timestamps(2:end-1)) ~= length(ttt - run.scan_start_ts)
                    disp(' aaa')
                    keyboard
                end
                if length(regs.state_timestamps) ~= length(regs.keystate_timestamps)
                    disp('w t f')
                    keyboard
                end
            end

        end

    end

    if ismember(regressors_collection, {'regressors', 'empa_regressors'})
        % to make sure the "visual" regressors from states (logged in core.py, in plays.states) line up with "EMPA" regressors (logged in EMPA.py, in regressors).
        % this is a challenge b/c we log ones in state timestamps and the others in keystate timestamps (logged separately in core.py), and b/c they have different lengths (b/c EMPA skips initial and last state)
        % so we have to manually adjust for that; in the end, we treat the keystates timestamps as ground truth
        assert(all(regs.state_timestamps > regs.keystate_timestamps)); % states are logged after keystates
        assert(mean(regs.state_timestamps - regs.keystate_timestamps) < 0.025); % ...but not too long after (~0.05 is the difference between two frames); notice if we're misaligned, lots of them will be off; maybe compare w/ mean(regs.state_timestamps(2:end) - regs.state_timestamps(1:end-1))

        regs.timestamps = regs.keystate_timestamps;

        X = [];

        for i = 1:numel(reg_fields)
            assert(size(regs.(reg_fields{i}),1) == length(regs.timestamps));
            X = [X, zscore(regs.(reg_fields{i}), 0, 1)];
        end
        for i = 1:numel(post_fields)
            assert(isequal(post_fields{i}, 'interaction_change_flag') || length(regs.(post_fields{i})) == length(regs.timestamps));
            X = [X, zscore(regs.(post_fields{i}), 0, 1)];
        end
    else

        X = [];
    end

    fields = [reg_fields, post_fields];

    if do_cache
        save(filename, 'regs', 'X', 'fields', '-v7.3');
    end
