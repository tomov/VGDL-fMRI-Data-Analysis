function show_figure(figure_name)


switch figure_name

    case 'plot_behavior'
        % plot_behavior.m
        %
        %load(fullfile(get_mat_dir(), 'plot_behavior.mat'));
        figure('pos', [370 777 462*0.6 435*0.6]) ;
        files = {'plot_behavior_1-11_dqn1200.mat', 'plot_behavior_12-32_dqn1200.mat'};

        % get subject scores for each agent (all groups lumped together)
        ys{1} = [];
        ys{2} = [];
        ys{3} = [];
        for group_idx = 1:length(files)
            load(fullfile(get_mat_dir(), 'show_figure', files{group_idx}), 'scores', 'wins', 'success', 'success_rates', 'plot_whats', 'agents', 'game_names');

            pw = 1; % expected payout
            plot_what = plot_whats{pw};

            for a = 1:length(ys)
                s = eval(plot_what);
                clear ss;
                for g = 1:length(game_names)
                    ss(:,g) = nanmean(s{g, a}, 2); % average across levels
                end
                s = mean(ss, 2); % and average across games too
                ys{a} = [ys{a}; s]; % lump together with other groups
            end
        end

        % Single plot for all groups
        %
        pw = 1; % expected payout
        plot_what = plot_whats{pw};

        agents = agents(1:3); % exclude random agent
        agent_center_offsets = [-0.25, 0, 0.25];

        %cmap = colormap(hsv(3));
        cmap = [0.4460 0.6740 0.1880;
                0 0.4470 0.7410;
                0.8500 0.3250 0.0980];

        % plots
        maxy = 0;
        for a = 1:length(agents)
            Violin(ys{a}, 0 + agent_center_offsets(a), 'ShowMean', true, 'Width', 0.1, 'ViolinColor', cmap(a, :));
            maxy = max(maxy, max(ys{a}));
        end

        % significance ***
        for a1 = 1:length(agents)
            for a2 = a1+1:length(agents)
                if all(isnan(ys{a1})) || all(isnan(ys{a2}))
                    continue
                end
                p = ranksum(ys{a1}, ys{a2}); % Mann Whitney U test across subjects
                if isnan(p), p = 1; end
                %maxy = maxy + 0.05;
                maxy = maxy + 0.8;
                y = maxy;
                x = mean([agent_center_offsets(a1) agent_center_offsets(a2)]);
                plot(0 + [agent_center_offsets(a1) agent_center_offsets(a2)], [y y], 'color', 'black');
                text(x, y + 0.2 + 0.2 * (p >= 0.05), significance(p), 'HorizontalAlignment', 'center', 'fontsize', 10);
            end
        end

        ylabel('Expected payout ($)');
        xlabel('Agent');
        set(gca, 'xtick', []);
        xlim([-0.38 0.38])
        ylim([0 14])

        % hacky custom legend
        h = zeros(length(agents), 1);
        for a = 1:length(agents)
            h(a) = plot(NaN, NaN, 'color', cmap(a,:));
        end
        %l = legend(h, {agents.name});
        if group_idx == length(files)
            l = legend(h, {'Human', 'EMPA', 'DDQN'});
            l.Position = [0.6651 0.6066 0.3234 0.1686];
        end

        title('Human and model behavior', 'interpreter', 'none');

        print('svg/figure2/plot_behavior.svg', '-dsvg');


    case 'plot_behavior_separate_DEPRECATED'
        % plot_behavior.m
        %
        %load(fullfile(get_mat_dir(), 'plot_behavior.mat'));
        figure('pos', [370 777 462*2*0.6 435*0.6]) ;
        files = {'plot_behavior_1-11_dqn1200.mat', 'plot_behavior_12-32_dqn1200.mat'};
        for group_idx = 1:length(files)
            load(fullfile(get_mat_dir(), 'show_figure', files{group_idx}));
            subplot(1, length(files), group_idx);

            pw = 1; % expected payout
            plot_what = plot_whats{pw};

            agents = agents(1:3); % exclude random agent
            agent_center_offsets = [-0.25, 0, 0.25];

            clear ys;
            maxy = 0;
            for a = 1:length(agents)
                s = eval(plot_what);
                clear ss;
                for g = 1:length(game_names)
                    ss(:,g) = nanmean(s{g, a}, 2); % average across levels
                end
                s = mean(ss, 2); % and average across games too
                Violin(s, 0 + agent_center_offsets(a), 'ShowMean', true, 'Width', 0.1, 'ViolinColor', cmap(a, :));
                ys{a} = s;
                maxy = max(maxy, max(ys{a}));
            end

            % significance ***
            for a1 = 1:length(agents)
                for a2 = a1+1:length(agents)
                    if all(isnan(ys{a1})) || all(isnan(ys{a2}))
                        continue
                    end
                    p = ranksum(ys{a1}, ys{a2}); % Mann Whitney U test across subjects
                    if isnan(p), p = 1; end
                    %maxy = maxy + 0.05;
                    maxy = maxy + 0.8;
                    y = maxy;
                    x = mean([agent_center_offsets(a1) agent_center_offsets(a2)]);
                    plot(0 + [agent_center_offsets(a1) agent_center_offsets(a2)], [y y], 'color', 'black');
                    text(x, y + 0.2 + 0.2 * (p >= 0.05), significance(p), 'HorizontalAlignment', 'center', 'fontsize', 10);
                end
            end

            ylabel('expected payout ($)');
            xlabel('agent');
            set(gca, 'xtick', []);
            xlim([-0.38 0.38])
            ylim([0 14])

            % hacky custom legend
            h = zeros(length(agents), 1);
            for a = 1:length(agents)
                h(a) = plot(NaN, NaN, 'color', cmap(a,:));
            end
            %l = legend(h, {agents.name});
            if group_idx == length(files)
                l = legend(h, {'Human', 'EMPA', 'DDQN'});
                l.Position = [0.3483 0.3806 0.1056 0.1149];
            end

            title(sprintf('Subjects %d..%d', min(subj_ids), max(subj_ids)), 'interpreter', 'none');
        end

        print('svg/figure2/plot_behavior_separate.svg', '-dsvg');

    otherwise
        assert(false, 'Invalid figure name');
end
