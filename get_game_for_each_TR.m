function [games, levels] = get_game_for_each_TR(subj_id)

    mongo_connect;

    nruns = 6;
    nTRs = 1698;
    initial_TRs = 7;
    n_games_per_run = 3;
    nTRs_per_game = (nTRs / nruns - initial_TRs) / n_games_per_run;

    EXPT = vgdl_expt;
    glmodel = 1;
    subj_id = 1;
    subj_games = {};

    % Get game order for subject

    for run_id = 1:nruns
        run = get_run(subj_id, run_id);

        [game_names, onsets, durs] = get_games(subj_id, run, conn);
        subj_games = [subj_games, game_names'];
    end

    % Expand games and levels to each TR

    games = {};
    levels = [];

    for r = 1:nruns
        games = [games; repmat({''}, [initial_TRs, 1])];
        levels = [levels; repmat(nan, [initial_TRs, 1])];
        
        p = partition_id_from_run_id(r);
        for g = (r - 1) * n_games_per_run + 1 : r * n_games_per_run
            games = [games; repmat(subj_games(g), [nTRs_per_game, 1])];
            levels = [levels; (ceil((1:nTRs_per_game) / nTRs_per_game * 3) + (p - 1) * 3)'];
        end
    end

    assert(length(games) == nTRs);
    assert(length(levels) == nTRs);
