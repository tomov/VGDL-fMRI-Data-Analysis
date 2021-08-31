clear all;
close all;

game_names = get_game_names_ordered(12);
subj_ids = 12:32;
%game_names = get_game_names_ordered(11);
%subj_ids = 1:11;

run_ids = 1:6;
levels = 1:9;

plot_what = 'success_rates'
assert(ismember(plot_what, {'scores', 'wins', 'success_rates'}));
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

% per level 

figure('pos', [64 421 2282 838]);

for g = 1:length(game_names)
    game_name = game_names{g};
    game_name
    
    x = [];
    y = [];

    scores{g} = nan(length(subj_ids), length(levels));
    wins{g} = nan(length(subj_ids), length(levels));
    success_rates{g} = nan(length(subj_ids), length(levels));
    for s = 1:length(subj_ids)
        subj_id = subj_ids(s);

        [instance_scores, instance_wins, instance_success_rates, instance_game_names, instance_levels] = get_instance_scores(subj_id, run_ids, true);
        which_instances = strcmp(instance_game_names, game_name);
        assert(sum(which_instances) >= 6);
        instance_scores = instance_scores(which_instances);
        instance_wins = instance_wins(which_instances);
        instance_success_rates = instance_success_rates(which_instances);
        instance_levels = instance_levels(which_instances);

        %assert(length(instance_scores) == length(levels)); -- not true if subject missed runs, e.g. subj 9
        assert(all(ismember(instance_levels, levels)));
        scores{g}(s, instance_levels) = instance_scores;
        wins{g}(s, instance_levels) = instance_wins;
        success_rates{g}(s, instance_levels) = instance_success_rates;
        
        x = [x instance_levels];
        y = [y eval(['instance_', plot_what])];
    end

    subplot(2, 3, g);
    hold on;

    %swarmchart(x, y);
    xs = cellstr(num2str(x'));
    violinplot(y, xs);
    xticks(levels);
    xlabel('level');
    ylabel(y_label);
    title(game_name, 'interpreter', 'none');
end


% per game

figure;

x = [];
y = [];

for g = 1:length(game_names)
    x = [x g * ones(1, length(subj_ids))];
    s = eval(plot_what);
    y = [y nanmean(s{g}, 2)'];  % average across levels for each subject
end

xs = cellstr(num2str(x'));
violinplot(y, xs);
xticks(1:length(game_names));
xticklabels(game_names);
set(gca,'TickLabelInterpreter','none');
xtickangle(30);
xlabel('game');
ylabel(y_label);
title(sprintf('Subjects %d..%d', min(subj_ids), max(subj_ids)), 'interpreter', 'none');
