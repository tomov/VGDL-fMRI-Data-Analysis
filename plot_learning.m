% copy of plot_behavior.m

clear all;
close all;

conn = mongo('holy7c22109.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54')

%game_names = get_game_names_ordered(11);
%subj_ids = 1:11;
game_names = get_game_names_ordered(12);
subj_ids = [12:32];
%subj_ids = [1:11];

run_ids = 1:6;
levels = 1:9;
frame_rate = 20; % Hz
level_duration = 60; %s
max_steps = frame_rate*level_duration*length(levels);
%
agents(1).name = 'Human';
agents(2).name = 'EMPA';
%agents(2).tag = 'attempt_1_states'; % 1..11
agents(2).tag = 'attempt_3_colors';  % 12..32
agents(3).name = 'DQN';
%agents(3).tag = '';   % 1..11
%%
%agents(4).name = 'Random';
%agents(4).tag = 'attempt_1_states';

%agents(3).name = 'EMPA';
%agents(3).tag = 'attempt_1_states';

%{
agents(1).name = 'EMPA';
agents(1).tag = 'attempt_3_colors';  % 12..32
agents(2).name = 'EMPA';
agents(2).tag = 'ablation_AGH3_attempt_1';  % 12..32
agents(3).name = 'EMPA';
agents(3).tag = 'ablation_IW_attempt_1';  % 12..32
agents(4).name = 'EMPA';
agents(4).tag = 'ablation_greedy_attempt_1';  % 12..32
agents(5).name = 'EMPA';
agents(5).tag = 'ablation_nodes_attempt_1';  % 12..32
agents(6).name = 'EMPA';
agents(6).tag = 'ablation_lessnodes_attempt_1';  % 12..32
%}


%% gather relevant data in a table

for a = 1:length(agents)
    %z = zeros(length(game_names), length(subj_ids), max_steps);
    z = zeros(length(game_names), length(subj_ids), max_steps);
    results(a).wins = z;
    results(a).cumwins = z;
end

human_steps = nan(length(game_names), length(subj_ids));

max_human_steps = 0;

for g = 1:length(game_names)
    game_name = game_names{g};
    game_name
        
    for a = 1:length(agents)
        agent_name = agents(a).name;
        agent_tag = agents(a).tag;
        agent_name

        for s = 1:length(subj_ids)
            subj_id = subj_ids(s);

            if strcmp(agent_name, 'Human')
                [play_scores, play_wins, play_actions, play_keypresses, play_durations, play_game_names, play_levels] = get_play_scores(conn, subj_id, run_ids, true);
                play_steps = play_actions;
                %play_steps = play_durations * frame_rate;
                %play_steps = play_keypresses;

            elseif strcmp(agent_name, 'DQN')
                %[play_scores, play_wins, play_steps, play_game_names, play_levels] = get_dqn_play_scores(game_name, subj_id, levels, true);
                %TODO!!!!!!!!!!!
                [play_scores, play_wins, play_steps, play_game_names, play_levels] = get_dqn_play_scores(game_name, 12, levels, true);
            else
                %[play_scores, play_wins, play_steps, play_game_names, play_levels] = get_agent_play_scores(conn, agent_name, subj_id, levels, agent_tag, true);
                %TODO!!!!!!!!!!!
                [play_scores, play_wins, play_steps, play_game_names, play_levels] = get_agent_play_scores(conn, agent_name, 12, levels, agent_tag, true);
            end

            % subselect game
            which_plays=strcmp(play_game_names,game_name); 

            play_wins=play_wins(which_plays);
            play_scores=play_scores(which_plays);
            play_steps=play_steps(which_plays);
            play_game_names=play_game_names(which_plays);
            play_levels=play_levels(which_plays);

            play_steps(play_steps == 0) = 1; % pad with an extra step when there were no steps in the episode
            cumsteps = cumsum(play_steps);

            if strcmp(agent_name, 'Human')
                human_steps(g,s) = max(cumsteps);
                max_human_steps = max(max_human_steps, max(cumsteps));
            end

            % truncate if there are too many steps
            play_wins(cumsteps > max_steps) = [];
            play_scores(cumsteps > max_steps) = [];
            play_steps(cumsteps > max_steps) = [];
            play_game_names(cumsteps > max_steps) = [];
            play_levels(cumsteps > max_steps) = [];
            cumsteps(cumsteps > max_steps) = [];

            % mark win events
            results(a).wins(g, s, cumsteps) = play_wins;
            results(a).cumwins(g, s, :) = cumsum(results(a).wins(g, s, :));

        end
    end
end


figure;

% each game separately yikes
for g = 1:length(game_names)
    game_name = game_names{g};
    game_name

    subplot(1,length(game_names),g);

    hold on;
    for a = 1:length(agents)
        data = squeeze(results(a).cumwins(g,:,:));
        m = mean(data,1);

        plot(m);
    end
    legend({agents.name});
    title(game_name,'interpreter','none');
end


figure;

% averaged across games
hold on;
for a = 1:length(agents)
    data = process_cumwins(results(a).cumwins,max_human_steps);

    m = mean(data,1);
    se = std(data,1)/sqrt(size(data,1));
    steps = 1:length(m);

    hp = plot(steps, m);
    hf = fill([steps flip(steps)], [m+se flip(m-se)], get(hp,'Color'));
    set(hf,'facealpha',0.3,'edgecolor','none');
end
legend({agents.name});
xlabel('steps');
ylabel('episodes won');
title('All games');


figure;

% averaged across games, separate subjects
hold on;
for a = 2:length(agents)
    data = process_cumwins(results(a).cumwins,max_human_steps);

    m = mean(data,1);
    se = std(data,1)/sqrt(size(data,1));
    steps = 1:length(m);

    plot(steps, m);
end

% plot humans

%data = process_cumwins(results(1).cumwins,max_human_steps);
load(fullfile(get_mat_dir(),'plot_learning_humans1-11_data.mat'), 'data');
data_1_11 = data;
load(fullfile(get_mat_dir(),'plot_learning_humans12-32_data.mat'), 'data');
data_12_32 = data;
clip = min(size(data_1_11,2),size(data_12_32,2));
data_1_11=data_1_11(:,1:clip);
data_12_32=data_12_32(:,1:clip);
data = [data_1_11; data_12_32];

plot(data', 'color', [0.5 0.5 0.5]);

legend('EMPA','DDQN','Human');
xlabel('steps');
ylabel('episodes won');
title('All games');
