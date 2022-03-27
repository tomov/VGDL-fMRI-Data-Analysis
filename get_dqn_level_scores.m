function [level_scores, level_wins, level_success, level_success_rates, game_names, actual_levels] = get_dqn_level_scores(game_name, subj_id, levels, do_cache)
    %agent_name = 'DQN';
    %subj_id = 13;
    %levels = [0, 1];

    % get score for each level (i.e. level) as max across all won plays 
    % this was the same way we determined to pay out for the fMRI study
    % same as get_level_scores.m except for artificial agents

    if ~exist('do_cache', 'var')
        do_cache = false;
    end
   
    %filename = fullfile(get_mat_dir(false), sprintf('get_dqn_level_scores_game=%s_subj=%d_levels=%s.mat', game_name, subj_id, sprintf('%d_', levels)));
    filename = fullfile(get_mat_dir(false), sprintf('get_dqn25m_level_scores_game=%s_subj=%d_levels=%s.mat', game_name, subj_id, sprintf('%d_', levels)));
    filename
    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    %filepath = fullfile('/n/home_fasse/mtomov13/RC_RL/reward_histories', sprintf('%s_reward_history_fmri_trial%d.csv', game_name, subj_id));
    %filepath = fullfile('/n/home_fasse/mtomov13/RC_RL/reward_histories_25M_eval_120000', sprintf('%s_reward_history_fmri_trial%d.csv', game_name, subj_id));
    %filepath = fullfile('/n/home_fasse/mtomov13/RC_RL/reward_histories_25M_eval_12000', sprintf('%s_reward_history_fmri_trial%d.csv', game_name, subj_id));
    filepath = fullfile('/n/home_fasse/mtomov13/RC_RL/reward_histories_25M_eval_1200', sprintf('%s_reward_history_fmri_trial%d.csv', game_name, subj_id));
    filepath
    T = readtable(filepath);

    level_scores = [];
    level_wins = [];
    level_success = [];
    level_success_rates = [];
    game_names = {};
    actual_levels = [];

    for l = 1:length(levels)
        level = levels(l);

        level_score = max(T.ep_reward((T.level == l - 1) & strcmp(T.win, 'True')));
        if isempty(level_score)
            level_score = 0;
        end
        %level_score = max(T.ep_reward(T.level == l - 1));
        level_win = sum(strcmp(T.win(T.level == l - 1), 'True'));
        nplays = sum(T.level == l - 1);

        level_scores = [level_scores, level_score];
        level_wins = [level_wins, level_win];
        level_success = [level_success, level_win > 0];
        level_success_rates = [level_success_rates, level_win / nplays];
        game_names = [game_names, game_name];
        actual_levels = [actual_levels, level];

    end
    if do_cache
        save(filename, 'level_scores', 'game_names', 'level_wins', 'level_success', 'level_success_rates', 'actual_levels', '-v7.3');
    end
