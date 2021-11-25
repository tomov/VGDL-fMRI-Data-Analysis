clear all;
close all;

conn = mongo('holy7c22101.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54', 'UserName', 'reader', 'Password', 'parolatamadafaqa')
%conn = mongo('holy7c22101.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54', 'UserName', 'root', 'Password', 'parolatabe')
%conn = mongo('holy7c22101.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54', 'UserName', 'root', 'Password', 'parolatabe', 'AuthMechanism', 'SCRAM_SHA_256')

%game_names = get_game_names_ordered(12);
%subj_ids = 12:32;
game_names = get_game_names_ordered(11);
subj_ids = 1:11;

run_ids = 1:6;
levels = 1:9;

agents(1).name = 'Human';
agents(2).name = 'Random';
agents(2).tag = 'attempt_1_states';
%agents(3).name = 'EMPA';
%agents(3).tag = 'attempt_1_states';

%plot_what = 'success_rates'
%plot_what = 'wins'
plot_what = 'scores'
plot_whats = {'scores', 'wins', 'success_rates'};
assert(ismember(plot_what, plot_whats));
y_label = get_plot_what_label(plot_what);

%% gather relevant data in a table

tbl_game_ix = [];
tbl_level = [];
tbl_agent_ix = [];
tbl_subj_id = [];
tbl_scores = [];
tbl_wins = [];
tbl_success_rates = [];

for g = 1:length(game_names)
    game_name = game_names{g};
    game_name
        
    for a = 1:length(agents)
        agent_name = agents(a).name;
        agent_tag = agents(a).tag;
        agent_name

        scores{g, a} = nan(length(subj_ids), length(levels));
        wins{g, a} = nan(length(subj_ids), length(levels));
        success_rates{g, a} = nan(length(subj_ids), length(levels));
        for s = 1:length(subj_ids)
            subj_id = subj_ids(s);

            if strcmp(agent_name, 'Human')
                [instance_scores, instance_wins, instance_success_rates, instance_game_names, instance_levels] = get_instance_scores(conn, subj_id, run_ids, true);
            else
                [instance_scores, instance_wins, instance_success_rates, instance_game_names, instance_levels] = get_agent_level_scores(conn, agent_name, subj_id, levels, agent_tag, true);
            end
            which_instances = strcmp(instance_game_names, game_name);
            %assert(sum(which_instances) >= 6);
            instance_scores = instance_scores(which_instances);
            instance_wins = instance_wins(which_instances);
            instance_success_rates = instance_success_rates(which_instances);
            instance_levels = instance_levels(which_instances);

            %assert(length(instance_scores) == length(levels)); -- not true if subject missed runs, e.g. subj 9
            assert(all(ismember(instance_levels, levels)));
            scores{g, a}(s, instance_levels) = instance_scores; % {game, agent}(subject, level) = average score
            wins{g, a}(s, instance_levels) = instance_wins;
            success_rates{g, a}(s, instance_levels) = instance_success_rates;

            tbl_game_ix = [tbl_game_ix, g * ones(1, length(instance_levels))];
            tbl_level = [tbl_level, instance_levels];
            tbl_agent_ix = [tbl_agent_ix, a * ones(1, length(instance_levels))];
            tbl_subj_id = [tbl_subj_id, subj_id * ones(1, length(instance_levels))];
            tbl_scores = [tbl_scores, instance_scores];
            tbl_wins = [tbl_wins, instance_wins];
            tbl_success_rates = [tbl_success_rates, instance_success_rates];
        end
    end

end

tbl = table(tbl_game_ix, tbl_level, tbl_agent_ix, tbl_subj_id, tbl_scores, tbl_wins, tbl_success_rates, 'VariableNames', {'game', 'level', 'agent', 'subj', 'score', 'win', 'success_rate'})


%% plot stuff
%

overall_width = 0.6; % width of violin plots for all agents
agent_width = 0.6 / length(agents); % width per agent, / 2 because Violin()
agent_center_offsets = - overall_width / 2 + agent_width * ((1:length(agents)) - 0.5); % center of each violin, for each agent

cmap = colormap(jet(length(agents)));


%% per level 

figure('pos', [64 421 2282 838]);

for g = 1:length(game_names)
    game_name = game_names{g};

    subplot(2, 3, g);
    hold on;

    for l = 1:length(levels)
        level = levels(l);

        clear ys;
        for a = 1:length(agents)
            s = eval(plot_what);
            s = s{g, a}(:, level);
            Violin(s, l + agent_center_offsets(a), 'ShowMean', true, 'Width', 0.1, 'ViolinColor', cmap(a, :));
            ys{a} = s;
        end

        % significance ***
        p = ranksum(ys{1}, ys{2}); % Mann Whitney U test across subjects
        if isnan(p), p = 1; end
        y = max([max(ys{1}) max(ys{2})]) + 0.5;
        plot(l + [agent_center_offsets(1) agent_center_offsets(2)], [y y], 'color', 'black');
        text(l, y + 0.3 + 0.3 * (p >= 0.05), significance(p), 'HorizontalAlignment', 'center', 'fontsize', 15);
    end

    % dividers
    ax = gca;
    for l = 0:length(levels)
        plot([l + 0.5, l + 0.5], ax.YLim, '--', 'color', [0.3 0.3 0.3]);
    end

    xticks(levels);
    xlabel('level');
    ylabel(y_label);
    title(game_name, 'interpreter', 'none');

    % hacky custom legend
    h = zeros(length(agents), 1);
    for a = 1:length(agents)
        h(a) = plot(NaN, NaN, 'color', cmap(a,:));
    end
    legend(h, {agents.name});

    hold off;
end



%% per game


figure;

hold on;

for g = 1:length(game_names)
    % average across levels for each subject
    clear ys;
    for a = 1:length(agents)
        s = eval(plot_what);
        s = nanmean(s{g, a}, 2); % average across levels
        Violin(s, g + agent_center_offsets(a), 'ShowMean', true, 'Width', 0.1, 'ViolinColor', cmap(a, :));
        ys{a} = s;
    end

    % significance ***
    p = ranksum(ys{1}, ys{2}); % Mann Whitney U test across subjects
    if isnan(p), p = 1; end
    y = max([max(ys{1}) max(ys{2})]) + 0.5;
    plot(g + [agent_center_offsets(1) agent_center_offsets(2)], [y y], 'color', 'black');
    text(g, y + 0.3 + 0.3 * (p >= 0.05), significance(p), 'HorizontalAlignment', 'center', 'fontsize', 15);
end

% dividers
ax = gca;
for g = 0:length(game_names)
    plot([g + 0.5, g + 0.5], ax.YLim, '--', 'color', [0.3 0.3 0.3]);
end

xticks(1:length(game_names));
xticklabels(game_names);
set(gca,'TickLabelInterpreter','none');
xtickangle(30);
xlabel('game');
ylabel(y_label);
title(sprintf('Subjects %d..%d', min(subj_ids), max(subj_ids)), 'interpreter', 'none');

% hacky custom legend
h = zeros(length(agents), 1);
for a = 1:length(agents)
    h(a) = plot(NaN, NaN, 'color', cmap(a,:));
end
legend(h, {agents.name});

hold off;


%% all games

figure('position', [420 917 1141 422]);

for pw = 1:length(plot_whats)
    plot_what = plot_whats{pw};

    subplot(1, length(plot_whats), pw);

    clear ys;
    for a = 1:length(agents)
        s = eval(plot_what);
        clear ss;
        for g = 1:length(game_names)
            ss(:,g) = nanmean(s{g, a}, 2); % average across levels
        end
        s = mean(ss, 2); % and average across games too
        Violin(s, 0 + agent_center_offsets(a), 'ShowMean', true, 'Width', 0.1, 'ViolinColor', cmap(a, :));
        ys{a} = s;
    end

    % significance ***
    p = ranksum(ys{1}, ys{2}); % Mann Whitney U test across subjects
    if isnan(p), p = 1; end
    y = max([max(ys{1}) max(ys{2})]) + 0.5;
    plot(0 + [agent_center_offsets(1) agent_center_offsets(2)], [y y], 'color', 'black');
    text(0, y + 0.1 + 0.3 * (p >= 0.05), significance(p), 'HorizontalAlignment', 'center', 'fontsize', 15);

    ylabel(get_plot_what_label(plot_what));
    set(gca, 'xtick', []);

    % hacky custom legend
    h = zeros(length(agents), 1);
    for a = 1:length(agents)
        h(a) = plot(NaN, NaN, 'color', cmap(a,:));
    end
    legend(h, {agents.name});

    title(sprintf('Subjects %d..%d', min(subj_ids), max(subj_ids)), 'interpreter', 'none');
end


function y_label = get_plot_what_label(plot_what)
    switch (plot_what)
        case 'scores'
            y_label = 'expected payout';
        case 'wins'
            y_label = '# wins / level';
        case 'success_rates'
            y_label = 'success rate';
        otherwise
            assert(false);
    end
end
