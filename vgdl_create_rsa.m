function rsa = vgdl_create_rsa(rsa_idx, subj_id, seed)

    save_output = true;

    % Create rsa structure, helper function for creating EXPT in
    % vgdl_expt.m 
    % copied from exploration_create_rsa.m
    %
    % USAGE: rsa = vgdl_create_rsa(model,subj_id)
    %
    % INPUTS:
    %   rsa_idx - positive integer indicating which RSA we're doing
    %   subj_id - integer specifying which subject is being analyzed
    %   seed (optional) - random seed that, if passed, is used to generate random permutation of labels, for permutation tests
    %
    % OUTPUTS:
    %   rsa - a structure with the following fields:
    %     .glmodel - which GLM to use to get the trial-by-trial betas; make sure to have a unique regressor for each trial, e.g. 'trial_onset_1', 'trial_onset_2', etc.
    %     .event - which within-trial event to use for neural activity; used to pick the right betas (needs to be substring of the regressor name), e.g. 'trial_onset'
    %     .mask - path to .nii file, or 3D binary vector of voxels to consider
    %     .radius - searchlight radius in voxels
    %     .which_betas - logical mask for which betas (trials) front the GLM to include (e.g. not timeouts)
    %     .model - struct array describing the models used for behavioral RDMs (see Kriegeskorte et al. 2008) with the fields:
    %         .name - model name
    %         .features - [nTrials x D] feature vector
    %         .runs - [nTrials x 1] run id vector
    %         .distance_measure - name (e.g. 'cosine') or function handler to be used as a distance measure for the RDMs (passed to pdist, see MATLAB documentation)
    %         .is_control - whether this is a control model (e.g. time)
    %
    % Momchil Tomov, May 2020

    fprintf('rsa %d, subj %d\n', rsa_idx, subj_id);

    [allSubjects, subjdirs, goodRuns, goodSubjs] = vgdl_getSubjectsDirsAndRuns();
    goodRun_ids = find(goodRuns{subj_id});
  
    filename = sprintf('mat/vgdl_create_rsa_%d_subj%d.mat', rsa_idx, subj_id)

    % hack to make it work on the cluster until they install MongoDB toolbox
    % just pre-generate the multi's locally, then load them on the cluster
    try
        conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54');
    catch e
        e
        fprintf('loading from %s\n', filename);
        load(filename);

        % special case b/c of cluster TODO better way...
        if exist('seed', 'var')
            % optionally shuffle game labels within each pair of runs, for permutation testing
            % this preserves all the essential structure of the experiment
            switch rsa_idx
                case 1
                    rng(seed + subj_id);
                    rsa.model(1).features = rsa.model(1).features([randperm(6) randperm(6)+6 randperm(6)+12]);
                otherwise
                    assert(false, 'invalid rsa_idx -- should be one of the above');
            end % end of switch statement
        end

        return
    end

    game_name_to_id = containers.Map({'vgfmri3_chase','vgfmri3_helper','vgfmri3_bait','vgfmri3_lemmings','vgfmri3_plaqueAttack','vgfmri3_zelda'}, 1:6);

    % RSAs
    %
    switch rsa_idx

        % basic with boxcars over blocks (i.e. GLM 1) 
        %
        case 1

            EXPT = vgdl_expt();

            features = [];
            runs = [];
            %goodRun_ids = goodRun_ids(randperm(length(goodRun_ids))); % TODO rm me -- more permutation testing, across runs
            for run_id = goodRun_ids
                % copied from vgdl_create_multi, GLM 1
                query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id); % in python we index runs from 0 (but not subjects) 

                run = find(conn, 'runs', 'query', query);
                assert(length(run) == 1);
                [game_names, onsets, durs] = get_games(subj_id, run, conn);

                %game_names = game_names(randperm(length(game_names))); % TOOD rm me -- permutation test, strong
                for i = 1:numel(game_names)
                    % one for each block
                    features = [features; game_name_to_id(game_names{i})];
                    runs = [runs; round(run_id/2)]; % lump adjacent runs together to exclude comparisons for trials within the same pair; otherwise we get biased in the negative direction, i.e. same games appear more different -- see below (same issue as when comparing within run, just for adjacent runs)
                end
                %features = [features; randperm(3)']; % TODO rm me -- more random test; basically, stuff nearby (same run or even neighboring runs) is correlated => we should exclude that stuff, or at least control for it; o/w we get lots of negative similarities (b/c all "same games" are far from each other, by design => will be different, compared to "different games", which tend to be closer to each other, on average)
            end
            %features = features(randperm(size(features, 1)), :); % TODO rm me -- permutation test, weak

            if exist('seed', 'var')
                % optionally shuffle game labels within each pair of runs, for permutation testing
                % this preserves all the essential structure of the experiment
                rng(seed + subj_id);
                features = features([randperm(6) randperm(6)+6 randperm(6)+12]);
            end

            rsa.event = 'vgfmri3_'; % just a prefix for the game name
            rsa.glmodel = 1;
            rsa.radius = 10 / 1.5;
            rsa.mask = 'masks/mask.nii';
            rsa.which_betas = logical(ones(size(features))); % all

            rsa.model(1).name = 'game';
            rsa.model(1).features = features;
            rsa.model(1).runs = runs;
            rsa.model(1).distance_measure = @(g1, g2) g1 ~= g2;
            rsa.model(1).is_control = false;


        % basic with beta series (i.e. deconvolved impulse regressors at every 2 s; see GLM 22)
        %
        case 2

            EXPT = vgdl_expt();

            features = [];
            runs = [];
            for run_id = goodRun_ids 
                % copied from vgdl_create_multi, GLM 22
                query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id); % in python we index runs from 0 (but not subjects) 

                run = find(conn, 'runs', 'query', query);
                assert(length(run) == 1);
                [game_names, onsets, durs] = get_games(subj_id, run, conn);

                for i = 1:numel(game_names)
                    for j = 1:length(onsets{i})
                        ons = [onsets{i}(j):EXPT.TR:onsets{i}(j)+durs{i}(j)];
                        for k = 1:length(ons)
                            % one for each TR / time point
                            features = [features; game_name_to_id(game_names{i})];
                            runs = [runs; round(run_id/2)]; % lump adjacent runs together to exclude comparisons for trials within the same pair; otherwise we get biased in the negative direction, i.e. same games appear more different -- see below (same issue as when comparing within run, just for adjacent runs)
                        end
                    end
                end
            end

            rsa.event = 'run_'; % just a prefix
            rsa.glmodel = 22;
            rsa.radius = 10 / 1.5;
            rsa.mask = 'masks/mask.nii';
            rsa.which_betas = logical(ones(size(features))); % all

            rsa.model(1).name = 'game';
            rsa.model(1).features = features;
            rsa.model(1).runs = runs;
            rsa.model(1).distance_measure = @(g1, g2) g1 ~= g2;
            rsa.model(1).is_control = false;

        otherwise
            assert(false, 'invalid rsa_idx -- should be one of the above');

    end % end of switch statement

    if save_output
        save(filename, 'rsa', '-v7.3'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
    end

    close(conn);

end
