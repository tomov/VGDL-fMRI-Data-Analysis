function [regs, X] = get_regressors(subj_id, run, conn)

    %{
    subj_id = 1; 
    run_id = 1; 
    query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id); 
    run = find(conn, 'runs', 'query', query);
    %}

    % helper function to get regressors for each frame in vgdl_create_multi
    % copied & modified from GLM 3
    % note run is a struct
    %

    assert(isequal(run.subj_id, num2str(subj_id)));

    % TODO interaction_change_flag when fixed
    %binreg_fields = {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag', 'newEffects_flag'}; % binary db.regressors => onsets only, durations irrelevant; have the option of having them as onsets only
    binreg_fields = {'theory_change_flag', 'sprite_change_flag', 'termination_change_flag', 'newEffects_flag'}; % binary db.regressors => onsets only, durations irrelevant; have the option of having them as onsets only
    reg_fields = [binreg_fields, {'likelihood', 'sum_lik', 'n_ts', 'num_effects', 'R_GG', 'R_GGs', 'R_SG', 'R_SGs'}]; % db.regressors 
    binpost_fields = {'interaction_change_flag'};
    % TODO d*_len when fixed
    %post_fields = [binpost_fields {'timestamps', 'S_len','I_len','T_len','Igen_len','Tnov_len','Ip_len','dS_len','dI_len','dT_len','dIgen_len','dTnov_len','dIp_len'}]; % db.plays_post 
    post_fields = [binpost_fields {'S_len','I_len','T_len','Igen_len','Tnov_len','Ip_len'}]; % db.plays_post 

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
    regs.timestamps = []; % these all come from EMPA, so the timestamps are keystate['ts'] timestamps, whereas those from the states have state['ts'] timestamps; TODO unify maybe

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
                assert(length(regressors) == 1);
                reg = regressors(1);

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
                regs.timestamps = [regs.timestamps; t - run.scan_start_ts]; % assuming all have the same ts


                plays_post = find(conn, 'plays_post', 'query', q);
                assert(length(plays_post) == 1);
                play_post = plays_post(1);

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
                    regs.(post_fields{i}) = [regs.(post_fields{i}); r];
                end
            end

        end

    end


    X = [];

    for i = 1:numel(reg_fields)
        assert(size(regs.(reg_fields{i}),1) == length(regs.timestamps));
        X = [X, zscore(regs.(reg_fields{i}), 0, 1)];
    end
    for i = 1:numel(post_fields)
        assert(isequal(post_fields{i}, 'interaction_change_flag') || length(regs.(post_fields{i})) == length(regs.timestamps));
        X = [X, zscore(regs.(post_fields{i}), 0, 1)];
    end
