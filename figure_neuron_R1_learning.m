clear all;

load(fullfile(get_mat_dir(false), 'neuron_R1_learning.mat'));

game_names = {'Chase','Helper','Bait','Lemmings',{'Plaque', 'Attack'}, {'Avoid', 'George'},'Zelda'};
game_ids = {1, 2, 3, 4, 5, 5, 6};
subj_ids = {1:32, 1:32, 1:32, 1:32, 1:11, 12:32, 1:32};
num_levels = 9;

%figure('pos', [99 96 1822 803]);
figure('pos', [99 181 1274 718]);
h = tiledlayout(length(game_names), num_levels, 'TileSpacing', 'none', 'Padding', 'none');

for g = 1:length(game_names)
    game_name = game_names{g};
    game_id = game_ids{g};
    subjs = subj_ids{g};

    for level = 1:num_levels
        tcf_smooth = nanmean(learning(game_id, level).tcf_smooth(subjs, :), 1);

        %subplot(length(game_names), num_levels, (game_id - 1) * num_levels + level);
        nexttile

        t = 1/frequency:1/frequency:level_duration;
        plot(t, tcf_smooth);

        %ylim([0 0.1]);
        %if level > 1
        %    set(gca, 'ytick', []);
        %end
        %if g < length(game_names)
        %    set(gca, 'xtick', [])
        %end
        if level == 1
            ylabel(game_name, 'fontweight','bold');
        end
        if g == length(game_names)
            xlabel('time (s)');
        end
        if g == 1
            title(sprintf('level %d', level));
        end
    end
end

sgtitle('Theory updates (smoothed)');

print('svg/neuron_revision/figure_neuron_R1_learning.svg', '-dsvg');
