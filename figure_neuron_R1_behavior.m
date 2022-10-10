% In response to neuron reviewer 1, comment about learning

function show_figure(figure_name)

switch figure_name
    
    case 'ablations'

        figure('pos', [370 777 300 400]) ;

        % baseline
        %
        filename = 'plot_behavior_12-32_dqn1200.mat';
        load(fullfile(get_mat_dir(), 'show_figure', filename), 'scores', 'wins', 'success', 'success_rates', 'plot_whats', 'agents', 'game_names');

        x_rows = [];
        s_rows = [];
        a = 1; % EMPA
        clear ss;
        for g = 1:length(game_names)
            ss(:,g) = nanmean(scores{g, a}, 2); % average across levels
        end
        s = mean(ss, 2); % and average across games too
        x_rows = [x_rows; repmat(1, length(s), 1)];
        s_rows = [s_rows; s];

        null = s;

        tbl = table(x_rows, s_rows, 'VariableNames', {'x', 's'});
        h = boxchart(tbl.x, tbl.s, 'BoxWidth', 1.0); 
        hold on;

        % ablations
        %
        filename = 'plot_behavior_12_ablations.mat';
        load(fullfile(get_mat_dir(), 'show_figure', filename), 'scores', 'wins', 'success', 'success_rates', 'plot_whats', 'agents', 'game_names');

        x = 2;
        xs = [];
        ys = [];
        for a = [2, 3, 5, 6]
            clear ss;
            for g = 1:length(game_names)
                ss(:,g) = nanmean(scores{g, a}, 2); % average across levels
            end
            s = mean(ss, 2); % and average across games too
            x_rows = [x_rows; x];
            xs = [xs; x];
            ys = [ys; s];
            s_rows = [s_rows; s];

            p = mean([null; s] < s);
            fprintf('%s -> p = %.4f\n', agents(a).tag, p);

            plot([x], [s], '.', 'MarkerSize', 30);
            x = x + 1;
        end
        
        hold off;
        xlim([0 6]);
        title('EMPA ablations');
        legend({'baseline', 'no intrinsic rewards', 'no pruning', '10x fewer nodes', '100x fewer nodes'});
        xticklabels({});
        ylabel('Expected payout ($)');
        
        print('svg/neuron_revision/ablations.svg', '-dsvg'); 

    otherwise
        assert(false, 'Invalid figure name');
end
