% optionally call after plot_behavior.m

%% per level 


for g = 1:length(game_names)
    game_name = game_names{g};

    subplot(2, 3, g);
    hold on;

    x = tbl.level(tbl.game == g & tbl.agent == 1);
    y = tbl.score(tbl.game == g & tbl.agent == 1);
    %swarmchart(x, y);
    xs = cellstr(num2str(x'));
    violinplot(y, xs);
    xticks(levels);
    xlabel('level');
    ylabel(y_label);
    title(game_name, 'interpreter', 'none');
end
%}


figure('pos', [64 421 2282 838]);

for g = 1:length(game_names)
    game_name = game_names{g};

    subplot(2, 3, g);
    hold on;

    for l = 1:length(levels)
        level = levels(l);

        %y = tbl.score(tbl.game == g & tbl.agent == 1 & tbl.level == level);
        %assert(length(y) == length(subj_ids));
        y = scores{g, 1}(:,level);
        Violin(y, level);
    end

    xticks(levels);
    xlabel('level');
    ylabel(y_label);
    title(game_name, 'interpreter', 'none');
end

%% per game

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

