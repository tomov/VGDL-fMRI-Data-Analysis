clear all;
close all;

game_names = get_game_names_ordered(12);
subj_ids = 12:32;
run_ids = 1:6;
levels = 1:9;

figure('pos', [64 421 2282 838]);

x = [];
y = [];

for g = 1:length(game_names)
    game_name = game_names{g};

    scores{g} = nan(length(subj_ids), length(levels));
    for s = 1:length(subj_ids)
        subj_id = subj_ids(s);

        [instance_scores, instance_wins, instant_success_rates, instance_game_names] = get_instance_scores(subj_id, run_ids, true);
        which_instances = strcmp(instance_game_names, game_name);
        instance_scores = instance_scores(which_instances);
        instance_wins = instance_wins(which_instances);
        instant_success_rates = instant_success_rates(which_instances);

        assert(length(instance_scores) == length(levels));
        x = [x levels];
        y = [y instance_scores];
    end

    subplot(2, 3, g);
    hold on;

    %swarmchart(x, y);
    xs = cellstr(num2str(x'));
    violinplot(y, xs);
    xticks(levels);
    xlabel('level');
    ylabel('expected payout');
    title(game_name, 'interpreter', 'none');
end
