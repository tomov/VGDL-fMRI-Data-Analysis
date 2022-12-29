% In response to neuron reviewer 1, comment about learning

function show_figure(figure_name)

switch figure_name
    
    case 'theory_update_timecourse'

        load(fullfile(get_mat_dir(false), 'neuron_R1_learning.mat'));

        game_names = {'Chase','Helper','Bait','Lemmings',{'Plaque', 'Attack'}, {'Avoid', 'George'},'Zelda'};
        game_ids = {1, 2, 3, 4, 5, 5, 6};
        subj_ids = {1:32, 1:32, 1:32, 1:32, 1:11, 12:32, 1:32};
        num_levels = 9;

        %figure('pos', [99 96 1822 803]);
        %figure('pos', [99 181 1274 718]);
        figure('pos', [99 201 1445 698]);
        h = tiledlayout(length(game_names), num_levels, 'TileSpacing', 'none', 'Padding', 'none');

        for g = 1:length(game_names)
            game_name = game_names{g};
            game_id = game_ids{g};
            subjs = subj_ids{g};

            for level = 1:num_levels
                tcf = learning(game_id, level).tcf_smooth(subjs, :);
                m = nanmean(tcf, 1);
                se = nanstd(tcf, 1) ./ sqrt(sum(~isnan(tcf), 1));

                %subplot(length(game_names), num_levels, (game_id - 1) * num_levels + level);
                nexttile

                t = 1/frequency:1/frequency:level_duration;

                hold on;
                plot(t, m);
                h = fill([t flip(t)], [m+se flip(m-se)], 'blue');
                set(h,'facealpha',0.3,'edgecolor','none');
                hold off;

                %ylim([0 0.2]);

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

        orient(gcf, 'landscape');
        print('svg/neuron_revision/figure_neuron_R1_learning_theory_update_timecourse.svg', '-dsvg');



    case 'theory_update_timecourse_merged'

        load(fullfile(get_mat_dir(false), 'neuron_R1_learning.mat'));

        game_names = {'Chase','Helper','Bait','Lemmings',{'Plaque', 'Attack'}, {'Avoid', 'George'},'Zelda'};
        game_ids = {1, 2, 3, 4, 5, 5, 6};
        subj_ids = {1:32, 1:32, 1:32, 1:32, 1:11, 12:32, 1:32};
        num_levels = 9;

        figure('pos', [99 201 1445 698]);
        h = tiledlayout(length(game_names), 1, 'TileSpacing', 'none', 'Padding', 'none');

        for g = 1:length(game_names)
            game_name = game_names{g};
            game_id = game_ids{g};
            subjs = subj_ids{g};

            tcfs = [];
            for level = 1:num_levels
                tcf_smooth = nanmean(learning(game_id, level).tcf_smooth(subjs, :), 1);
                tcfs(level, :) = tcf_smooth;
            end

            nexttile

            t = 1/frequency:1/frequency:level_duration;
            plot(t, tcfs);
            xlabel('time (s)');
            ylabel(game_name, 'fontweight','bold');
            legend({'level 1', 'level 2', 'level 3', ...
                    'level 4', 'level 5', 'level 6', ...
                    'level 7', 'level 8', 'level 9'});
        end

        sgtitle('Theory updates (smoothed)');

        orient(gcf, 'landscape');
        print('svg/neuron_revision/figure_neuron_R1_learning_theory_update_timecourse_merged.svg', '-dsvg');



    case 'theory_update_timecourse_2'

        load(fullfile(get_mat_dir(false), 'neuron_R1_learning_2.mat'));
        %load(fullfile(get_mat_dir(false), 'neuron_R1_learning_2_sigma=400.000.mat'));

        runs_per_petition=2;

        tcf = [learning.tcf_smooth]';
        m = nanmean(tcf, 1);
        se = nanstd(tcf, 1) ./ sqrt(sum(~isnan(tcf), 1));
        t = 1/frequency : 1/frequency : run_duration * num_runs;

        plot(t, m);
        %h = fill([t flip(t)], [m+se flip(m-se)], 'blue');
        %set(h,'facealpha',0.3,'edgecolor','none');
        %plot(tcf(11, :)');

        snaohtu
        keyboard

        % run/data partition annotations
        for r=1:num_runs
            t_prev=(r-1)*run_duration;
            t_next=r*run_duration;
            if (r==1)
                plot([t_prev t_prev], [-0.2 1.0], '--', 'color', [0.5 0.5 0.5]);
            end
            plot([t_next t_next], [-0.2 1.0], '--', 'color', [0.5 0.5 0.5]);
            text(mean([t_prev t_next]), 0.8, sprintf('run %d', r), 'fontsize', 10, 'HorizontalAlignment', 'center');
        end
        for p=1:num_runs/runs_per_partition
            t_prev=(p-1)*run_duration*jruns_per_petition;
            t_next=p*run_duration*runs_per_petition;
            if (p==1)
                plot([t_prev t_prev], [-0.2 1.3], '--', 'color', [0.5 0.5 0.5]);
            end
            plot([t_next t_next], [-0.2 1.3], '--', 'color', [0.5 0.5 0.5]);
            text(mean([t_prev t_next]), 1.1, sprintf('partition %d', p), 'fontsize', 10, 'HorizontalAlignment', 'center');
        end

        ylabel('theory updates (a.u.)');
        xlabel('time (s)');

        title('Theory updates (smoothed) for entire session');

        orient(gcf, 'landscape');
        print('svg/neuron_revision/figure_neuron_R1_learning_theory_update_timecourse_2.svg', '-dsvg');



    case 'theory_update_timecourse_3'

        load(fullfile(get_mat_dir(false), 'neuron_R1_learning.mat'));
        %load(fullfile(get_mat_dir(false), 'neuron_R1_learning_sigma=400.mat'));

        game_names = {'Chase','Helper','Bait','Lemmings','Plaque Attack', 'Avoid George','Zelda'};
        game_ids = {1, 2, 3, 4, 5, 5, 6};
        subj_ids = {1:32, 1:32, 1:32, 1:32, 1:11, 12:32, 1:32};
        num_levels = 9;
        levels_per_partition = 3;
        num_partitions = 3;

        %figure('pos', [99 96 1822 803]);
        %figure('pos', [99 181 1274 718]);
        figure('pos', [99 201 1445 398]);
        hh = tiledlayout(2, 4, 'TileSpacing', 'none', 'Padding', 'none');

        % data wrangling
        for g = 1:length(game_names)
            game_name = game_names{g};
            game_id = game_ids{g};
            subjs = subj_ids{g};

            for subj = subjs
                tcfs_by_gs{g,subj} = [];
                for level = 1:num_levels
                    tcfs_by_gs{g,subj} = [tcfs_by_gs{g,subj}; learning(game_id, level).tcf_smooth(subj, :)'];
                end
            end
            tcf_by_g_by_s{g} = [tcfs_by_gs{g,:}]';
        end
        for subj = subjs
            tcf_by_s(subj,:) = nanmean([tcfs_by_gs{:,subj}]', 1);
        end

        for g = 1:length(game_names) + 1

            nexttile;
            hold on;

            % theory update histograms
            if g <= length(game_names)
                game_name = game_names{g};
                game_id = game_ids{g};
                subjs = subj_ids{g};

                data = tcf_by_g_by_s{g}(subjs,:);
                m = nanmean(data, 1);
                %se = nanstd(tcf_by_g_by_s{g}, 1) ./ sqrt(size(tcf_by_g_by_s{g}, 1));
                se = nanstd(data, 1) ./ sqrt(sum(~isnan(data), 1));
                t = 1/frequency:1/frequency:level_duration*num_levels;

                plot(t,m);
                %plot(t,tcf_by_g_by_s{g}');
                h = fill([t flip(t)], [m+se flip(m-se)], 'blue');
                set(h,'facealpha',0.3,'edgecolor','none');

                xlim([0 level_duration*num_levels]);
                xlabel('time (s)');
                ylabel('theory updates (a.u.)');
                title(game_name);

            else
                % all games (special case)

                m = nanmean(tcf_by_s, 1);
                se = nanstd(tcf_by_s, 1) ./ sqrt(sum(~isnan(tcf_by_s), 1));
                t = 1/frequency:1/frequency:level_duration*num_levels;

                plot(t,m);
                h = fill([t flip(t)], [m+se flip(m-se)], 'blue');
                set(h,'facealpha',0.3,'edgecolor','none');

                xlim([0 level_duration*num_levels]);
                xlabel('time (s)');
                ylabel('theory updates (a.u.)');
                title('All games');
            end

            % level/data partition annotations
            yl=ylim;
            for level=1:num_levels
                t_prev=(level-1)*level_duration;
                t_next=(level)*level_duration;
                if(level==1)
                    plot([t_prev t_prev], [yl(1) yl(2)*0.8], '--','color',[0.5 0.5 0.5]); 
                end
                plot([t_next t_next], [yl(1) yl(2)*0.8], '--','color',[0.5 0.5 0.5]); 
                text(mean([t_prev t_next]), yl(2)*0.7, sprintf('l%d', level), 'fontsize', 10, 'HorizontalAlignment', 'center');
            end
            for part=1:num_partitions
                t_prev=(part-1)*level_duration*levels_per_partition;
                t_next=(part)*level_duration*levels_per_partition;
                if(level==1)
                    plot([t_prev t_prev], [yl(1) yl(2)*1.0], '-','color',[0.5 0.5 0.5]); 
                end
                plot([t_next t_next], [yl(1) yl(2)*1.0], '-','color',[0.5 0.5 0.5]); 
                text(mean([t_prev t_next]), yl(2)*0.9, sprintf('partition %d', part), 'fontsize', 10, 'HorizontalAlignment', 'center');
            end


            hold off;
        end

        sgtitle('Theory updates (smoothed)');

        orient(gcf, 'landscape');
        print('svg/neuron_revision/figure_neuron_R1_learning_theory_update_timecourse_3.svg', '-dsvg');



    case 'avn_timecourse'

        load(fullfile(get_mat_dir(false), 'neuron_R1_learning_avn_sigma=20.mat'));
        %load(fullfile(get_mat_dir(false), 'neuron_R1_learning_sigma=400.mat'));

        game_names = {'Chase','Helper','Bait','Lemmings','Plaque Attack', 'Avoid George','Zelda'};
        game_ids = {1, 2, 3, 4, 5, 5, 6};
        subj_ids = {1:32, 1:32, 1:32, 1:32, 1:11, 12:32, 1:32};

        %figure('pos', [99 96 1822 803]);
        %figure('pos', [99 181 1274 718]);
        figure('pos', [99 201 1445 398]);
        hh = tiledlayout(2, 4, 'TileSpacing', 'none', 'Padding', 'none');

        for g = 1:length(game_names) + 1

            nexttile;
            hold on;

            % theory update histograms
            if g <= length(game_names)
                game_name = game_names{g};
                game_id = game_ids{g};
                subjs = subj_ids{g};

                for v=1:numel(valences)

                    data = learning_avn_smooth(g).(valences{v}); 
                    m = nanmean(data, 1);
                    %se = nanstd(tcf_by_g_by_s{g}, 1) ./ sqrt(size(tcf_by_g_by_s{g}, 1));
                    se = nanstd(data, 1) ./ sqrt(sum(~isnan(data), 1));
                    %t = 1/frequency:1/frequency:level_duration*num_levels;

                    plot(m);
                    %plot(t,m);
                    %plot(t,tcf_by_g_by_s{g}');
                    %h = fill([t flip(t)], [m+se flip(m-se)], 'blue');
                    %set(h,'facealpha',0.3,'edgecolor','none');
                end

                legend(valences);
                %xlim([0 level_duration*num_levels]);
                xlabel('time (s)');
                ylabel('theory updates (a.u.)');
                title(game_name);

            else
                % all games (special case)

                %m = nanmean(tcf_by_s, 1);
                %se = nanstd(tcf_by_s, 1) ./ sqrt(sum(~isnan(tcf_by_s), 1));
                %t = 1/frequency:1/frequency:level_duration*num_levels;

                %plot(t,m);
                %h = fill([t flip(t)], [m+se flip(m-se)], 'blue');
                %set(h,'facealpha',0.3,'edgecolor','none');

                %xlim([0 level_duration*num_levels]);
                %xlabel('time (s)');
                %ylabel('theory updates (a.u.)');
                %title('All games');
            end


            hold off;
        end

        sgtitle('sntoheu (smoothed)');

        orient(gcf, 'landscape');
        print('svg/neuron_revision/figure_neuron_R1_learning_theory_update_timecourse_3.svg', '-dsvg');





    case 'GLM_102'

        ccnl_view(vgdl_expt(), 102, 'theory_change_flag');

    case 'GLM_196_first_partition'

        ccnl_view(vgdl_expt(), 196, 'theory_change_flag_part1');

    case 'GLM_196_second_partition'

        ccnl_view(vgdl_expt(), 196, 'theory_change_flag_part2');

    case 'GLM_196_third_partition'

        ccnl_view(vgdl_expt(), 196, 'theory_change_flag_part3');

    case 'GLM_196_contrast'

        ccnl_view(vgdl_expt(), 196, 'theory_change_flag_part1 - theory_change_flag_part3');

    case 'plot_gp_CV_EMPA_partitions_1_2'
        % plot_gp_CV.m

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=12.mat');
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_EMPA_partitions_2_3'
        % plot_gp_CV.m

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=23.mat');
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);


    case 'plot_gp_CV_EMPA_partitions_contrast'
        % plot_gp_CV.m

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=12.mat');
        zs_12 = zs;

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=23.mat');
        zs_23 = zs;

        [h,p,ci,stats] = ttest(zs_12, zs_23);
        ts = stats.tstat;
        tmap(mask) = ts;

        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    otherwise
        assert(false, 'Invalid figure name');
end
