%function multi = vgdl_create_multi(glmodel, subj_id, run_id, save_output)
glmodel = 1;
subj_id = 181;
run_id = 1;

clear multi;

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
    run_ids = find(goodRuns{find(subj_id == allSubjects)});
    run_id = run_ids(run_id);
    fprintf('run_id %d \n', run_id);

    conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')

    query = sprintf('{"subj_id": "%d"}', subj_id)

    subj = find(conn, 'subjects', 'query', query)
    assert(length(subj) > 0, 'Subject not found -- wrong subj_id?');
    assert(length(subj) <= 1, 'Too many subjects -- duplicate entries? Yikes!');

    query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id - 1) % in python we index runs from 0 (but not subjects) 

    run = find(conn, 'runs', 'query', query)
    assert(length(run) == 1);
    
    plays = find(conn, 'plays', 'query', query)

    regressors = find(conn, 'regressors', 'query', query)
    assert(length(plays) == length(regressors));

    % GLMs
    %
    switch glmodel

        % condition = game, boxcars over blocks
        % goal is to look for systematic differences across games for things like agency, etc
        %
        case 1 

            idx = 0;

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
                    en = instance.end_time - run.scan_start_ts;

                    onsets = [onsets, st];
                    durs = [durs, en - st];
                end

                idx = idx + 1;
                multi.names{idx} = game_name;
                multi.onsets{idx} = onsets;
                multi.durations{idx} = durs;
            end

            % TODO add events -- visuals & key down / key ups

        case 2
            multi.names{1} = 'test';

        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');
    end % end of switch statement

    if save_output
        save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
    end


