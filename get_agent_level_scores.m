function [level_scores, level_wins, level_success_rates, game_names] = get_agent_level_scores(agent_name, subj_id, levels, do_cache)
    %agent_name = 'DQN';
    %subj_id = 13;
    %levels = [0, 1];

    % get score for each level (i.e. level) as max across all won plays 
    % this was the same way we determined to pay out for the fMRI study
    % same as get_level_scores.m except for artificial agents

    if ~exist('do_cache', 'var')
        do_cache = false;
    end
   
    filename = sprintf('mat/get_agent_level_scores_agent=%s_subj=%d_levels=%s.mat', agent_name, subj_id, sprintf('%d_', levels));
    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    level_scores = [];
    level_wins = [];
    level_success_rates = [];
    game_names = {};

    conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')
    for l = 1:length(levels)
        level = levels(l);

        q = sprintf('{"subj_id": "%d", "agent_name": "%s", "results.level": %d}', subj_id, agent_name, level);
        res = find(conn, 'sim_results', 'query', q);
        for r = 1:length(res)
            % multiple levels from the same game
            % only consider the level that we care about
            
            level_score = 0;
            level_win = 0;
            nplays = 0;
            game_name = res(r).game_name;

            for rr = 1:length(res(r).results)
                sim_play = res(r).results(rr); 

                if sim_play.level ~= level
                    continue
                end

                if sim_play.win
                    level_score = max(level_score, sim_play.score);
                    level_win = level_win + 1;
                end
                nplays = nplays + 1;
            end

            level_scores = [level_scores, level_score];
            level_wins = [level_wins, level_win];
            level_success_rates = [level_success_rates, level_win / nplays];
            game_names = [game_names, game_name];
        end
    end

    if do_cache
        save(filename, 'level_scores', 'game_names', 'level_wins', 'level_success_rates', '-v7.3');
    end
