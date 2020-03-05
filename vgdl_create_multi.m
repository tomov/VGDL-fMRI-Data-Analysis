function multi = vgdl_create_multi(glmodel, subj, run, save_output)

    % Create multi structure, helper function for creating EXPT in
    % vgdl_expt.m
    % copied from exploration_create_multi.m from https://github.com/tomov/Exploration-Data-Analysis
    %
    % USAGE: multi = vgdl_create_multi(model,subj,run)
    %
    % INPUTS:
    %   glmodel - positive integer indicating general linear model
    %   subj - integer specifying which subject is being analyzed
    %   run - integer specifying the run
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

    if nargin < 4 || isempty(save_output)
        save_output = false;
    end

    fprintf('glm %d, subj %d, run %d\n', glmodel, subj, run);

    [allSubjects, subjdirs, goodRuns, goodSubjs] = vgdl_getSubjectsDirsAndRuns();

    SPM_run = run; % save the SPM run (SPM doesn't see bad runs)
  
    % skip bad runs
    runs = find(goodRuns{subj});
    run = runs(run);
    fprintf('run %d \n', run);
   

    % GLMs
    %
    switch glmodel

        case 1 
            multi.names{1} = 'test';

        case 2
            multi.names{1} = 'test';

        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');
    end % end of switch statement

    if save_output
        save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
    end
end
