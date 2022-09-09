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

        orient(gcf, 'landscape');
        print('svg/neuron_revision/figure_neuron_R1_learning.svg', '-dsvg');
        print('pdf/neuron_revision/figure_neuron_R1_learning.pdf', '-dpdf', '-bestfit');

    case 'GLM_102'

        ccnl_view(vgdl_expt(), 102, 'theory_change_flag');

    case 'GLM_196_first_partition'

        ccnl_view(vgdl_expt(), 196, 'theory_change_flag_part1');

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
