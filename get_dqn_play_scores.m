function [play_scores, play_wins, play_steps, game_names, actual_levels] = get_dqn_play_scores(game_name, subj_id, levels, do_cache)
    % copy of get_dqn_level_scores.m but for plays/episodes

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
    filename = fullfile(get_mat_dir(false), sprintf('get_dqn25m_play_scores_game=%s_subj=%d_levels=%s.mat', game_name, subj_id, sprintf('%d_', levels)));
    filename
    if do_cache
        if exist(filename, 'file')
            load(filename);
            return
        end
    end

    filepath = fullfile('/n/home_fasse/mtomov13/RC_RL/reward_histories_25M_eval_1200', sprintf('%s_reward_history_fmri_trial%d.csv', game_name, subj_id));
    filepath
    T = readtable(filepath);

    play_scores = [];
    play_wins = [];
    play_steps = [];
    game_names = {};
    actual_levels = [];

    last_steps = 0;
    for p = 1:size(T,1)
        if ismember(T.level(p) + 1, levels)
            play_scores = [play_scores, T.ep_reward(p)];
            play_wins = [play_wins, strcmp(T.win(p), 'True')];
            play_steps = [play_steps, T.steps(p)-last_steps];
            game_names = [game_names, game_name];
            actual_levels = [actual_levels, T.level(p) + 1];
        end
        last_steps = T.steps(p);
    end

    if do_cache
        save(filename, 'play_scores','play_wins','play_steps','game_names','actual_levels', '-v7.3');
    end
