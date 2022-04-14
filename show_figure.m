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

        print('pdf/plot_behavior.pdf', '-dpdf');
        print('svg/plot_behavior.svg', '-dsvg');


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

        print('pdf/plot_behavior_separate.pdf', '-dpdf');
        print('svg/plot_behavior_separate.svg', '-dsvg');

    %
    %
    %% Figure 3: EMPA theory vs. DDQN layer representations
    %
    %

    case 'plot_gp_CV_EMPA'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1.mat')); % this is it !!!!!!!!!!!!!!!!!!!!ontroller five
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_DDQN'
        % plot_gp_CV.m

        %load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_all_nsamples=100_project=1_norm=1_fast=1.mat')); % !!!
        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN25M_all_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1.mat')); % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_PCA'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=PCA_all_nsamples=100_project=1_norm=1_fast=1.mat')); % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_VAE'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=VAE_e1k__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1.mat')); % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_VAE_e10k'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=VAE__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1.mat')); % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_game'
        % plot_gp_CV.m
        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=game__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1.mat'));
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_EMPA_noproject'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1.mat')); % 
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_DDQN_noproject'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN25M_all_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1.mat')); % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_PCA_noproject'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=PCA__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1.mat')); % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA'
        % plot_gp_CV_rois.m

        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_vae.mat');
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_vae_repro.mat');
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_25M_e1k.mat');
        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA_25M_e1k.mat');
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.001_atlas=AAL2_GP_EMPA_25M_e1k.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 521 1045 218]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});
        %h = plot_gp_CV_rois_helper(fs(:,ix,:), 'ranksum', 'median', regressor_names(ix), roi_names, [], [], 1, [1:2]);
        %h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 2, [1:1]);
        %h = plot_gp_CV_rois_helper(zs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names, [], [], 1, [1:2]);
        h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 8, [1:1], 1.2);
        title('Model comparison by ROI');
        ylabel('Fraction significant voxels');

        % Prettify it
        yscale = 7.2;
        xscale = 1.25;
        text(xscale * 3.5, yscale * 0.075, 'Frontal/Motor', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 6.5 xscale * 6.5], yscale * [0 0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 8.5, yscale * 0.075, 'Dorsal/Parietal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 10.5 xscale * 10.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 12.6, yscale *0.075, 'Ventral/Temporal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 14.5 xscale * 14.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 15.6, yscale * 0.075, 'Early visual', 'fontsize', 10, 'HorizontalAlignment', 'center');
        ylim([0 yscale *0.08]);
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        l = legend({'EMPA', 'DDQN', 'PCA', 'VAE'});
        l.Position = [0.1371 0.6119 0.0842 0.2821];

        orient(gcf, 'landscape');
        print('pdf/plot_gp_CV_rois_fraction_AAL2_GP_EMPA.pdf', '-dpdf');
        print('svg/plot_gp_CV_rois_fraction_AAL2_GP_EMPA.svg', '-dsvg');


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped'
        % plot_gp_CV_rois.m

        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped.mat');
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped_vae.mat');
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped_vae_repro.mat');
        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA_grouped_25M_e1k.mat');
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped_25M_e1k.mat');
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.001_atlas=AAL2_GP_EMPA_grouped_25M_e1k.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 522 537 417]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});
        %h = plot_gp_CV_rois_helper(fs(:,ix,:), 'ranksum', 'median', regressor_names(ix), roi_names, [], [], 1, [1:2]);
        h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 2, [1:1]);
        %h = plot_gp_CV_rois_helper(fs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names, [], [], 1, [1:2]);
        title('Model comparison by ROI group');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'});
        l = legend({'EMPA', 'DDQN', 'PCA', 'VAE'});

        print('pdf/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped.pdf', '-dpdf');


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped_EMPA_components'
        % plot_gp_CV_rois.m

        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped.mat');
        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA_grouped_25M_e1k.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 522 537 417]);
        ix = ismember(regressor_names, {'sprite', 'interaction', 'termination'});
        cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap, 2); %colormap(winter(3)));
        title('EMPA components');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'});
        l = legend({'objects', 'relations', 'goals'});
        l.Position = [0.1799 0.7430 0.1922 0.1205];
        ylim([0 0.25]);

        print('pdf/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped_EMPA_components.pdf', '-dpdf');

    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped_DDQN_layers'
        % plot_gp_CV_rois.m

        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped.mat');
        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA_grouped_25M_e1k.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 522 537 417]);
        ix = ismember(regressor_names, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
        cmap = [0.9 0.8 0.7 0.6 0.5]' * [0.8500    0.3250    0.0980] + [0.1 0.2 0.3 0.4 0.5]' * [1 1 1];
        plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap, 2); %colormap(autumn(5)));
        title('DDQN layers');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'});
        l = legend({'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
        l.Position = [0.1626 0.6843 0.1732 0.1942];

        print('pdf/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped_DDQN_layers.pdf', '-dpdf');


    case 'plot_ridge_CV_EMPA'
        % plot_ridge_CV.m

        load(fullfile(get_mat_dir(), 'agg_ridge_CV_us=1_glm=1_model=EMPA_theory_subsample=0_project=1.mat'));
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    %
    %
    %% Figure 4: EMPA theory updating
    %
    %


    case 'GLM_102'

        ccnl_view(vgdl_expt(), 102, 'theory_change_flag');

    case 'plot_confirmatory_betas_for_masks_AAL2_GLM_102'
        % plot_confirmatory_betas_for_masks.m

        load(fullfile(get_mat_dir(false), 'confirmatory_betas_for_masks_atlas=AAL2_GLM_102_cglm=102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-.mat')); % !!!!!!!!!!!!!!!!!!!!!
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        betas = betas(ROI_ix);

        subjs = 1:1:32;

        figure('position', [97 451 2190 888]);

        %confirmatory_regressors = confirmatory_regressors(end-2:end);
        ps = [];

        for m = 1:length(mask_filenames)
        %for m = 1:2
            beta = betas{m}(subjs, :);
            %beta = beta(:,end-2:end);
            beta

            [sem, me] = wse(beta);
            [h,p,ci,stats] = ttest(beta);

            subplot(4,4,m);
            %subplot(1,2,m);
            hold on;
            bar(me);
            h = errorbar(me, sem, 'o', 'MarkerSize', 2);
            h.CapSize = 1;

            ax = gca;
            for j = 1:size(beta, 2)
                if p(j) <= 0.05
                    %text(j, ax.YLim(2) - 0.1 + 0.1 * rand, significance(p(j)), 'fontsize', 3, 'HorizontalAlignment', 'center'); 
                    if me(j) < 0
                        y = me(j) - sem(j) - 0.1;
                    else
                        y = me(j) + sem(j) + 0.1;
                    end
                    text(j, y, significance(p(j)), 'fontsize', 3, 'HorizontalAlignment', 'center'); 
                end
            end

            ylabel('beta coefficient');
            ax.TickLabelInterpreter = 'none';
            xticklabels(confirmatory_regressors);
            set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',4)
            xtickangle(30);
            xticks(1:length(confirmatory_regressors));
            title(regions{m}, 'interpreter', 'none', 'fontsize', 10);
            ylim(1.1 * ylim);

            %p_corr = 1 - (1 - p(1)) ^ length(mask_filenames);
            %text(10, mean([0 ax.YLim(2)]), sprintf('p_{corr.} = %.2e', p_corr), 'fontsize', 17);
            ps(m) = p(1);
        end

        ps_corr = bonferroni(ps)';
        table(regions, mask_name', ps', ps_corr)

        regions(ps_corr < 0.05)

        ROI_ix = find(ps_corr <= 0.05)

        h = gcf;
        set(h,'PaperOrientation','landscape');
        print('pdf/plot_confirmatory_betas_for_masks_AAL2_GLM_102.pdf', '-dpdf', '-fillpage');


    case 'plot_confirmatory_betas_for_masks_AAL2_GLM_102_components'
        % plot_confirmatory_betas_for_masks.m

        load(fullfile(get_mat_dir(false), 'confirmatory_betas_for_masks_atlas=AAL2_GLM_102_cglm=103-104-105-.mat')); % !!!!!!!!!
        ROI_ix = [      1      2      8     11     12     13     14     15   16]; 

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        betas = betas(ROI_ix);

        subjs = 1:1:32;

        figure('position', [97 451 2190 888]);

        %confirmatory_regressors = confirmatory_regressors(end-2:end);
        ps = [];

        for m = 1:length(mask_filenames)
        %for m = 1:2
            beta = betas{m}(subjs, :);
            %beta = beta(:,end-2:end);
            beta

            [sem, me] = wse(beta);
            [h,p,ci,stats] = ttest(beta);

            subplot(4,4,m);
            %subplot(1,2,m);
            hold on;
            bar(me);
            errorbar(me, sem, 'o', 'MarkerSize', 1);

            ax = gca;
            for j = 1:size(beta, 2)
                if p(j) <= 0.05
                    text(j, ax.YLim(2) - 0.1, significance(p(j)), 'fontsize', 10, 'HorizontalAlignment', 'center'); 
                end
            end

            ylabel('beta coefficient');
            ax.TickLabelInterpreter = 'none';
            xticklabels(confirmatory_regressors);
            set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',8)
            xtickangle(30);
            xticks(1:length(confirmatory_regressors));
            title(regions{m}, 'interpreter', 'none');

            %p_corr = 1 - (1 - p(1)) ^ length(mask_filenames);
            %text(10, mean([0 ax.YLim(2)]), sprintf('p_{corr.} = %.2e', p_corr), 'fontsize', 17);
            ps(m) = p(1);
        end

        ps_corr = bonferroni(ps)';
        table(regions, mask_name', ps', ps_corr)

        regions(ps_corr < 0.05)

        ROI_ix = find(ps_corr <= 0.05)

    case 'plot_PETH_AAL2_GLM_102_BOLD'
        % plot_PETHs.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GLM_102_BOLD.mat')); % !!!!!!!!!!!!
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [      1      2      8     11     12     13     14     15   16]; 

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        figure('pos', [64 348 2110*0.5 911*0.5]);

        % optionally plot theory change flag only
        %fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        fields(find(strcmp(fields, 'termination_change_flag'))) = [];
        fields(find(strcmp(fields, 'block_start'))) = [];
        fields(find(strcmp(fields, 'block_end'))) = [];
        fields(find(strcmp(fields, 'instance_start'))) = [];
        fields(find(strcmp(fields, 'instance_end'))) = [];

        subjs = 1:1:32;

        cmap = colormap(jet(length(fields)));
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:length(mask_filenames)
            disp(mask_name{m});

            subplot(3, 4, m + (floor((m - 1) / 3)));
            hold on;

            for i = 1:length(fields)
                field = fields{i};
                disp(field)

                D = activations(m).(field)(subjs,:); % subj x TRs PETH's
                [sem, me] = wse(D);
                %me = nanmean(activations(m).(field), 1);
                %sem = nanstd(activations(m).(field), 1) / sqrt(size(activations(m).(field), 1)); % TODO wse
                [h,p,ci,stats] = ttest(D);

                %errorbar(dTRs, m, se);
                hh(i) = plot(t, me, 'color', cmap(i,:));
                h = fill([t flip(t)], [me + sem flip(me - sem)], cmap(i,:));
                set(h, 'facealpha', 0.3, 'edgecolor', 'none');

                ax = gca;
                xlim([t(1) - 1, t(end) + 1]);
                if ismember(field, regs_fields)
                    ix = find(strcmp(field, regs_fields)) - 1;
                    j_init = find(PETH_dTRs > 0); % ignore baselines
                    for j = j_init:length(t)
                        if p(j) <= 0.05
                            h = text(t(j) + ix * 0.1, ax.YLim(2) - 0.02 - ix * 0.02, significance(p(j)), 'color', cmap(i,:), 'fontsize', 5, 'HorizontalAlignment', 'center');
                            set(h,'Rotation',90);
                        end
                    end
                end

                plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
                plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            end

            if m == length(mask_filenames)
                %legend(hh, fields, 'interpreter', 'none');
                l = legend(hh, {'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
                l.Position = [0.7399 0.7786 0.0948 0.1570];
                %l.Position = [0.0196 0.5936 0.0946 0.1573];
            end
            if exist('what', 'var') && strcmp(what, 'GP')
                ylabel('\Delta z');
            else
                ylabel('\Delta BOLD');
            end
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
        end

        orient(gcf, 'landscape');
        %print('pdf/plot_PETH_AAL2_GLM_102_BOLD.pdf', '-dpdf', '-fillpage');
        print('pdf/plot_PETH_AAL2_GLM_102_BOLD.pdf', '-dpdf');


    case 'plot_PETHs_bars_AAL2_GLM_102_BOLD'
        % plot_PETHs_bars.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GLM_102_BOLD.mat')); % !!!!!!!!!!!!
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [      1      2      8     11     12     13     14     15   16]; 

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        % optionally plot theory change flag only
        %fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        fields(find(strcmp(fields, 'termination_change_flag'))) = [];
        fields(find(strcmp(fields, 'block_start'))) = [];
        fields(find(strcmp(fields, 'block_end'))) = [];
        fields(find(strcmp(fields, 'instance_start'))) = [];
        fields(find(strcmp(fields, 'instance_end'))) = [];

        subjs = 1:1:32;

        nROIs = length(mask_filenames);
        nregressors = length(fields);
        nsubjects = length(subjs);

        as = nan(nROIs,nregressors,nsubjects);

        cmap = colormap(jet(length(fields)));
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:nROIs
            disp(mask_name{m});

            for i = 1:nregressors
                field = fields{i};
                disp(field)

                D = activations(m).(field)(subjs,:); % subj x TRs PETH's
                as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
            end
        end

        % Piggyback off of plot_gp_CV_rois.m
        figure('pos', [49 329 2143 610]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 9, 1:1);
        if exist('what', 'var') && strcmp(what, 'GP')
            title('Average Fisher z-transformed Pearson correlation change');
            ylabel('\Delta z');
        else
            title('Average BOLD change in ROIs');
            ylabel('\Delta BOLD');
        end

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.75, 'Frontal/Motor', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(3.05, 0.75, 'Dorsal/Parietal', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([3.5 3.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(5.5, 0.75, 'Ventral/Temporal', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([7.5 7.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(8.5, 0.75, 'Early visual', 'fontsize', 9, 'HorizontalAlignment', 'center');
        %legend(fields(ix), 'interpreter', 'none');
        l = legend({'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
        l.Position = [0.6830 0.6710 0.0742 0.2344];
        l.FontSize = 6;

        % fix the significant stars
        H = findobj(gcf);
        H(15).Position=[9.3429 0.1004 0];

        orient(gcf, 'landscape');
        print('pdf/plot_PETHs_bars_AAL2_GLM_102_BOLD.pdf', '-dpdf', '-bestfit');



    case 'plot_PETH_components_AAL2_GLM_102_BOLD'
        % plot_PETHs.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GLM_102_BOLD.mat')); % !!!!!!!!!!!!
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [      1      2      8     11     12     13     14     15   16]; 

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        figure('pos', [64 348 2110*0.5 911*0.5]);

        % optionally plot theory change flag only
        fields = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};
        %fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        %%fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        %%fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        %%fields(find(strcmp(fields, 'termination_change_flag'))) = [];
        %fields(find(strcmp(fields, 'block_start'))) = [];
        %fields(find(strcmp(fields, 'block_end'))) = [];
        %fields(find(strcmp(fields, 'instance_start'))) = [];
        %fields(find(strcmp(fields, 'instance_end'))) = [];

        subjs = 1:1:32;

        %cmap = colormap(jet(length(fields)));
        cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:length(mask_filenames)
            disp(mask_name{m});

            subplot(3, 4, m + (floor((m - 1) / 3)));
            hold on;

            for i = 1:length(fields)
                field = fields{i};
                disp(field)

                D = activations(m).(field)(subjs,:); % subj x TRs PETH's
                [sem, me] = wse(D);
                %me = nanmean(activations(m).(field), 1);
                %sem = nanstd(activations(m).(field), 1) / sqrt(size(activations(m).(field), 1)); % TODO wse
                [h,p,ci,stats] = ttest(D);

                %errorbar(dTRs, m, se);
                hh(i) = plot(t, me, 'color', cmap(i,:));
                h = fill([t flip(t)], [me + sem flip(me - sem)], cmap(i,:));
                set(h, 'facealpha', 0.3, 'edgecolor', 'none');

                ax = gca;
                xlim([t(1) - 1, t(end) + 1]);
                if ismember(field, regs_fields)
                    ix = find(strcmp(field, regs_fields)) - 1;
                    j_init = find(PETH_dTRs > 0); % ignore baselines
                    for j = j_init:length(t)
                        if p(j) <= 0.05
                            h = text(t(j) + ix * 0.5, ax.YLim(2) - 0.02 - ix * 0.02, significance(p(j)), 'color', cmap(i,:), 'fontsize', 7, 'HorizontalAlignment', 'center');
                            set(h,'Rotation',90);
                        end
                    end
                end

                plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
                plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            end

            if m == length(mask_filenames)
                %legend(hh, fields, 'interpreter', 'none');
                %l = legend(hh, {'object update', 'relation update', 'goal update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
                l = legend(hh, {'object update', 'relation update', 'goal update'}, 'interpreter', 'none');
                l.Position = [0.7399 0.7786 0.0948 0.1570];
            end
            if exist('what', 'var') && strcmp(what, 'GP')
                ylabel('\Delta z');
            else
                ylabel('\Delta BOLD');
            end
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
        end

        orient(gcf, 'landscape');
        %print('pdf/plot_PETH_components_AAL2_GLM_102_BOLD.pdf', '-dpdf', '-fillpage');
        print('pdf/plot_PETH_components_AAL2_GLM_102_BOLD.pdf', '-dpdf');


    case 'plot_PETHs_bars_components_AAL2_GLM_102_BOLD'
        % plot_PETHs_bars.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GLM_102_BOLD.mat')); % !!!!!!!!!!!!
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [      1      2      8     11     12     13     14     15   16]; 

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        % optionally plot theory change flag only
        fields = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};
        %fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        %%fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        %%fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        %%fields(find(strcmp(fields, 'termination_change_flag'))) = [];
        %fields(find(strcmp(fields, 'block_start'))) = [];
        %fields(find(strcmp(fields, 'block_end'))) = [];
        %fields(find(strcmp(fields, 'instance_start'))) = [];
        %fields(find(strcmp(fields, 'instance_end'))) = [];

        subjs = 1:1:32;

        nROIs = length(mask_filenames);
        nregressors = length(fields);
        nsubjects = length(subjs);

        as = nan(nROIs,nregressors,nsubjects);

        %cmap = colormap(jet(length(fields)));
        cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:nROIs
            disp(mask_name{m});

            for i = 1:nregressors
                field = fields{i};
                disp(field)

                D = activations(m).(field)(subjs,:); % subj x TRs PETH's
                as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
            end
        end

        % Piggyback off of plot_gp_CV_rois.m
        figure('pos', [49 329 943 410]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 8, 1:3);
        if exist('what', 'var') && strcmp(what, 'GP')
            title('Average Fisher z-transformed Pearson correlation change');
            ylabel('\Delta z');
        else
            title('Average BOLD change in ROIs');
            ylabel('\Delta BOLD');
        end

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.50, 'Frontal/Motor', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(3.05, 0.50, 'Dorsal/Parietal', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([3.5 3.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(5.5, 0.50, 'Ventral/Temporal', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([7.5 7.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(8.5, 0.50, 'Early visual', 'fontsize', 9, 'HorizontalAlignment', 'center');
        %legend(fields(ix), 'interpreter', 'none');
        %l = legend({'object update', 'relation update', 'goal update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
        l = legend({'object update', 'relation update', 'goal update'}, 'interpreter', 'none');
        l.Position = [0.6830 0.6710 0.0742 0.0944];
        l.FontSize = 6;
        ylim([-0.2 0.55]);


        orient(gcf, 'landscape');
        %print('pdf/plot_PETHs_bars_components_AAL2_GLM_102_BOLD.pdf', '-dpdf', '-bestfit');
        print('pdf/plot_PETHs_bars_components_AAL2_GLM_102_BOLD.pdf', '-dpdf');



    case 'plot_glm_bic_bms_AAL2_GLM_102_multiplex_with_controls'
        % plot_glm_bic_bms.m

        load(fullfile(get_mat_dir(false), 'glm_bic_bms_atlas=AAL2_GLM_102_multiplex_with_controls.mat'));
        glm_ix = [1 2 3 4 5];
        ROI_ix = [      1      2      8     11     12     13    ]; 
        mask_filenames = mask_filenames(ROI_ix);
        regions = regions(ROI_ix);
        bics = bics(ROI_ix);

        subjs = 1:1:32;

        nROIs = length(mask_filenames);
        nGLMs = length(glms);
        nsubjects = length(subjs);

        bs = nan(nROIs,nGLMs,nsubjects);

        clear pxps;
        clear bors;

        for m = 1:length(mask_filenames)
            %bic = bics{m}(subjs, [2 3 4 8]);
            %bic = bics{m}(subjs, [1 2 3 4 8]);
            [~, mask_name{m}, ~] = fileparts(mask_filenames{m});
            bic = bics{m}(subjs, glm_ix);
            bs(m,glm_ix,subjs) = bic';

            lme = -0.5 * bic;
            [alpha, exp_r, xp, pxp, bor] = bms(lme);
            pxps(m,:) = pxp;
            bors(m,:) = bor;
        end

        % Show table
        subjs
        glm_names(glm_ix)
        table(regions, pxps, bors)

        bs = bs - bs(:,1,:); % subtract first model

        % Piggyback off of plot_gp_CV_rois.m
        figure('pos', [49 568 1252 371]);

        %cmap = [1 0.8 0.6 0.4 0.2]' * [0    0.4470    0.7410];
        %cmap = [ ...
        % 0         0    1.0000; ...
        % 0    0.3333    1.0000; ...
        % 0    0.6667    1.0000];
        cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        cmap = [cmap; 0.7294    0.3333    0.8275];
        glm_names = {'theory updates', 'object updates', 'relation updates', 'goal updates', 'object, relation, goal updates'};
        ix = [2,3,4,5];
        h = plot_gp_CV_rois_helper(bs(:,ix,:), 'ttest', 'mean', glm_names(ix), regions, 0, cmap, 11, []);
        ylabel('\Delta BIC relative to theory updates');
        title('GLM comparison');

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.85 * 50000, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 1.8 * 50000], '--', 'color', [0.5 0.5 0.5]);
        text(3.05, 0.85 * 50000, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([3.5 3.5], [0 1.8 * 50000], '--', 'color', [0.5 0.5 0.5]);
        text(5, 0.85 * 50000, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        ylim([-40000 50000]);
        l = legend(glm_names(ix),  'interpreter', 'none');
        l.Position = [0.2619 0.2435 0.1815 0.1775];
        l.FontSize = 8;

        % minor adjustments
        H = findobj(gcf);
        H(12).Position(2) = H(12).Position(2) - 3000;
        H(16).Position(2) = H(16).Position(2) - 3000;
        H(19).Position(2) = H(19).Position(2) - 3000;
        H(26).Position(2) = H(26).Position(2) - 3000;

        orient(gcf, 'landscape');
        print('pdf/plot_glm_bic_bms_AAL2_GLM_102_multiplex_with_controls.pdf', '-dpdf', '-bestfit');

    %
    %
    %% Figure 5
    %
    %


    case 'contrast_overlap_GP_EMPA_GLM_102'
        % contrast_overlap.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1.mat'));
        tmap_filename = bspmview_save_map(EXPT, tmap);

        EXPTs = {vgdl_expt(), tmap_filename}; % Take advantage of hack that allows us to plot any custom nmap
        glmodels = [102 , -1];
        contrasts = {'theory_change_flag', ''};

        % group-level settings
        p = 0.001;
        alpha = 0.05;
        Dis = 20;
        if ~exist('extent', 'var')
            extent = []; % use default cluster size
        end
        if ~exist('Num', 'var')
            Num = 1; % # peak voxels per cluster; default in bspmview is 3
        end
        direct = '+';
        df = 31;
        clusterFWEcorrect = true;

        for i = 1:length(glmodels)
            EXPT = EXPTs{i};
            glmodel = glmodels(i);
            contrast = contrasts{i};
            extent = [];

            [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent, df);

            if ~exist('cumulative_mask', 'var')
                cumulative_mask = C > 0;
            else
                cumulative_mask = cumulative_mask + (C>0);
            end
        end

        bspmview_wrapper(vgdl_expt(), cumulative_mask);


    case 'plot_PETH_AAL2_GP_EMPA_GLM_102_GP'
        % plot_PETHs.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP.mat')); % !!!!!!!!!!!!
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_CV_no_baseline.mat')); 
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        figure('pos', [64 460 1296*0.5 799*0.5]);

        % optionally plot theory change flag only
        %fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        fields(find(strcmp(fields, 'termination_change_flag'))) = [];
        fields(find(strcmp(fields, 'block_start'))) = [];
        fields(find(strcmp(fields, 'block_end'))) = [];
        fields(find(strcmp(fields, 'instance_start'))) = [];
        fields(find(strcmp(fields, 'instance_end'))) = [];

        subjs = 1:1:32;

        cmap = colormap(jet(length(fields)));
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:length(mask_filenames)
            disp(mask_name{m});

            subplot(3, 3, m);
            hold on;

            for i = 1:length(fields)
                field = fields{i};
                disp(field)

                D = activations(m).(field)(subjs,:); % subj x TRs PETH's
                [sem, me] = wse(D);
                %me = nanmean(activations(m).(field), 1);
                %sem = nanstd(activations(m).(field), 1) / sqrt(size(activations(m).(field), 1)); % TODO wse
                [h,p,ci,stats] = ttest(D);

                %errorbar(dTRs, m, se);
                hh(i) = plot(t, me, 'color', cmap(i,:));
                h = fill([t flip(t)], [me + sem flip(me - sem)], cmap(i,:));
                set(h, 'facealpha', 0.3, 'edgecolor', 'none');

                ax = gca;
                xlim([t(1) - 1, t(end) + 1]);
                if ismember(field, regs_fields)
                    ix = find(strcmp(field, regs_fields)) - 1;
                    j_init = find(PETH_dTRs > 0); % ignore baselines
                    for j = j_init:length(t)
                        if p(j) <= 0.05
                            if mean(me) < 0
                                y = ax.YLim(1) + 0.01 - ix * 0.01;
                            else
                                y = ax.YLim(2) - 0.01 - ix * 0.01;
                            end
                            h = text(t(j) + ix * 0.1, y, significance(p(j)), 'color', cmap(i,:), 'fontsize', 7, 'HorizontalAlignment', 'center');
                            set(h,'Rotation',90);
                        end
                    end
                end

                plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
                plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            end

            if m == 1
                %l = legend(hh, fields, 'interpreter', 'none');
                l = legend(hh, {'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
                l.Position = [0.4023 0.1654 0.1289 0.1314];
            end
            %ylabel('z');
            ylabel('\Delta z');
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
        end

        orient(gcf, 'landscape');
        %print('pdf/plot_PETH_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf', '-bestfit');
        print('pdf/plot_PETH_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf');
 

    case 'plot_PETH_components_AAL2_GP_EMPA_GLM_102_GP'
        % plot_PETHs.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_sprite.mat')); % !!
        sprite_activations = activations;
        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_interaction.mat')); % !!
        interaction_activations = activations;
        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_termination.mat')); % !!
        termination_activations = activations;
        clear activations;
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP.mat')); % !!!!!!!!!!!!
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_CV_no_baseline.mat')); 
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        sprite_activations = sprite_activations(ROI_ix);
        interaction_activations = interaction_activations(ROI_ix);
        termination_activations = termination_activations(ROI_ix);

        figure('pos', [64 460 1296*0.5 799*0.5]);

        % only look at the components
        fields = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

        subjs = 1:1:32;

        %cmap = colormap(jet(length(fields)));
        %cmap = colormap(jet(9));
        cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:length(mask_filenames)
            disp(mask_name{m});

            subplot(3, 3, m);
            hold on;

            for i = 1:length(fields)
                field = fields{i};
                disp(field)

                % hardcoded
                if  i == 1
                    D = sprite_activations(m).(field)(subjs,:); % subj x TRs PETH's
                elseif i == 2
                    D = interaction_activations(m).(field)(subjs,:); % subj x TRs PETH's
                elseif i == 3
                    D = termination_activations(m).(field)(subjs,:); % subj x TRs PETH's
                else
                    assert(false, 'we are matching component update events with the corresponding component activations; unclear which activations to show for non-component update events');
                end
                [sem, me] = wse(D);
                %me = nanmean(activations(m).(field), 1);
                %sem = nanstd(activations(m).(field), 1) / sqrt(size(activations(m).(field), 1)); % TODO wse
                [h,p,ci,stats] = ttest(D);
                %errorbar(dTRs, m, se);
                hh(i) = plot(t, me, 'color', cmap(i,:));
                h = fill([t flip(t)], [me + sem flip(me - sem)], cmap(i,:));
                set(h, 'facealpha', 0.3, 'edgecolor', 'none');

                ax = gca;
                xlim([t(1) - 1, t(end) + 1]);
                if ismember(field, regs_fields)
                    ix = find(strcmp(field, regs_fields)) - 1;
                    j_init = find(PETH_dTRs > 0); % ignore baselines
                    for j = j_init:length(t)
                        if p(j) <= 0.05
                            h = text(t(j) + ix * 0.5, ax.YLim(2) + 0.005 + ix * 0.002, significance(p(j)), 'color', cmap(i,:), 'fontsize', 7, 'HorizontalAlignment', 'center');
                            set(h,'Rotation',90);
                        end
                    end
                end

                plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
                plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            end

            if m == 1
                %legend(hh, fields, 'interpreter', 'none');
                l = legend(hh, {'object update', 'relation update', 'goal update'}, 'interpreter', 'none');
                l.Position = [0.4162 0.1941 0.1073 0.0588];
            end
            ylabel('\Delta z');
            %ylabel('z');
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
            ll = ylim;
            ylim([ll(1) ll(2)*1.2]);
        end

        orient(gcf, 'landscape');
        %print('pdf/plot_PETH_components_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf', '-bestfit');
        print('pdf/plot_PETH_components_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf');


    case 'plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP'
        % plot_PETHs_bars.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP.mat')); % !!!!!!!!!!!!
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_CV_no_baseline.mat')); 
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        % optionally plot theory change flag only
        %fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        fields(find(strcmp(fields, 'termination_change_flag'))) = [];
        fields(find(strcmp(fields, 'block_start'))) = [];
        fields(find(strcmp(fields, 'block_end'))) = [];
        fields(find(strcmp(fields, 'instance_start'))) = [];
        fields(find(strcmp(fields, 'instance_end'))) = [];

        subjs = 1:1:32;

        nROIs = length(mask_filenames);
        nregressors = length(fields);
        nsubjects = length(subjs);

        as = nan(nROIs,nregressors,nsubjects);

        cmap = colormap(jet(length(fields)));
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:nROIs
            disp(mask_name{m});

            for i = 1:nregressors
                field = fields{i};
                disp(field)

                D = activations(m).(field)(subjs,:); % subj x TRs PETH's
                as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
                %as(m,i,:) = mean(D(:, PETH_dTRs > 5), 2); % average across time, ignoring baseline
            end
        end

        % Piggyback off of plot_gp_CV_rois.m
        figure('pos', [49 329 2143*0.5 610*0.5]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 2, 1, 4);
        title('Average Fisher z-transformed Pearson correlation change');
        ylabel('\Delta z');
        %title('Average Fisher z-transformed Pearson correlation in ROIs');
        %ylabel('z');

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.15, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.3], '--', 'color', [0.5 0.5 0.5]);
        text(3.5, 0.15, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([4.5 4.5], [0 0.3], '--', 'color', [0.5 0.5 0.5]);
        text(6, 0.15, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        %l = legend(fields(ix), 'interpreter', 'none');
        l = legend({'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
        l.Position = [0.5597 0.6938 0.0807 0.1836];
        l.FontSize = 7;
        ylim([-0.05 0.2]);

        orient(gcf, 'landscape');
        %print('pdf/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf', '-bestfit');
        print('pdf/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf');



    case 'plot_PETH_bars_AAL2_GP_EMPA_GLM_102_grouped_GP'
        % plot_PETHs_bars.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_grouped_GP.mat')); 
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_grouped_GP_CV_no_baseline.mat')); 
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        % optionally plot theory change flag only
        %fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        fields(find(strcmp(fields, 'termination_change_flag'))) = [];
        fields(find(strcmp(fields, 'block_start'))) = [];
        fields(find(strcmp(fields, 'block_end'))) = [];
        fields(find(strcmp(fields, 'instance_start'))) = [];
        fields(find(strcmp(fields, 'instance_end'))) = [];

        subjs = 1:1:32;

        nROIs = length(mask_filenames);
        nregressors = length(fields);
        nsubjects = length(subjs);

        as = nan(nROIs,nregressors,nsubjects);

        cmap = colormap(jet(length(fields)));
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:nROIs
            disp(mask_name{m});

            for i = 1:nregressors
                field = fields{i};
                disp(field)

                D = activations(m).(field)(subjs,:); % subj x TRs PETH's
                as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
                %as(m,i,:) = mean(D(:, PETH_dTRs > 5), 2); % average across time, ignoring baseline
            end
        end

        % Piggyback off of plot_gp_CV_rois.m
        figure('pos', [47 487 456 248]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 2, 1, 6);
        title('Average Fisher z-transformed Pearson correlation change');
        ylabel('\Delta z');
        %title('Average Fisher z-transformed Pearson correlation in ROIs');
        %ylabel('z');

        % Prettyfy it 
        %legend(fields(ix), 'interpreter', 'none');
        l = legend({'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
        l.Position = [0.3663 0.6003 0.2637 0.1600];
        l.FontSize = 6;
        xticklabels({'Frontal/Motor (IFG)', 'Dorsal/Parietal', 'Ventral/Temporal'});
        ylim([-0.02 0.125]);

         
        orient(gcf, 'landscape');
        print('pdf/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_grouped_GP.pdf', '-dpdf');
         

    case 'plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP_components'
        % plot_PETHs_bars.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_sprite.mat')); % !!
        sprite_activations = activations;
        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_interaction.mat')); % !!
        interaction_activations = activations;
        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_termination.mat')); % !!
        termination_activations = activations;
        clear activations;
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP.mat')); % !!!!!!!!!!!!
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_CV_no_baseline.mat')); 
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        sprite_activations = sprite_activations(ROI_ix);
        interaction_activations = interaction_activations(ROI_ix);
        termination_activations = termination_activations(ROI_ix);

        fields = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

        subjs = 1:1:32;

        nROIs = length(mask_filenames);
        nregressors = length(fields);
        nsubjects = length(subjs);

        as = nan(nROIs,nregressors,nsubjects);

        %cmap = colormap(jet(9));
        cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:nROIs
            disp(mask_name{m});

            for i = 1:nregressors
                field = fields{i};
                disp(field)

                % hardcoded
                if  i == 1
                    D = sprite_activations(m).(field)(subjs,:); % subj x TRs PETH's
                elseif i == 2
                    D = interaction_activations(m).(field)(subjs,:); % subj x TRs PETH's
                elseif i == 3
                    D = termination_activations(m).(field)(subjs,:); % subj x TRs PETH's
                else
                    assert(false, 'we are matching component update events with the corresponding component activations; unclear which activations to show for non-component update events');
                end

                as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
                %as(m,i,:) = mean(D(:, PETH_dTRs > 5), 2); % average across time, ignoring baseline
            end
        end

        % Piggyback off of plot_gp_CV_rois.m
        figure('pos', [750 464 1457 471]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 3, 1:3, 4);
        title('Average Fisher z-transformed Pearson correlation change');
        ylabel('\Delta z');
        %title('Average Fisher z-transformed Pearson correlation in ROIs');
        %ylabel('z');

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.25, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.35], '--', 'color', [0.5 0.5 0.5]);
        text(3.5, 0.25, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([4.5 4.5], [0 0.35], '--', 'color', [0.5 0.5 0.5]);
        text(6, 0.25, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        %legend(fields(ix), 'interpreter', 'none');
        ylim([-0.05 0.29]);
        l = legend({'object update', 'relation update', 'goal update'}, 'interpreter', 'none');
        l.Position = [0.5404 0.5824 0.0988 0.1062];

        orient(gcf, 'landscape');
        print('pdf/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP_components.pdf', '-dpdf', '-bestfit');



    case 'plot_PETH_bars_AAL2_GP_EMPA_GLM_102_grouped_GP_components'
        % plot_PETHs_bars.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_grouped_GP_sprite.mat')); % !!
        sprite_activations = activations;
        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_grouped_GP_interaction.mat')); % !!
        interaction_activations = activations;
        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_grouped_GP_termination.mat')); % !!
        termination_activations = activations;
        clear activations;
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_grouped_GP.mat')); 
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        sprite_activations = sprite_activations(ROI_ix);
        interaction_activations = interaction_activations(ROI_ix);
        termination_activations = termination_activations(ROI_ix);

        % optionally plot theory change flag only
        fields = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

        subjs = 1:1:32;

        nROIs = length(mask_filenames);
        nregressors = length(fields);
        nsubjects = length(subjs);

        as = nan(nROIs,nregressors,nsubjects);

        cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:nROIs
            disp(mask_name{m});

            for i = 1:nregressors
                field = fields{i};
                disp(field)

                % hardcoded
                if  i == 1
                    D = sprite_activations(m).(field)(subjs,:); % subj x TRs PETH's
                elseif i == 2
                    D = interaction_activations(m).(field)(subjs,:); % subj x TRs PETH's
                elseif i == 3
                    D = termination_activations(m).(field)(subjs,:); % subj x TRs PETH's
                else
                    assert(false, 'we are matching component update events with the corresponding component activations; unclear which activations to show for non-component update events');
                end

                as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
                %as(m,i,:) = mean(D(:, PETH_dTRs > 5), 2); % average across time, ignoring baseline
            end
        end

        % Piggyback off of plot_gp_CV_rois.m
        figure('pos', [1489 437 703 502]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 1, 1:3, 5);
        title('Average Fisher z-transformed Pearson correlation change');
        ylabel('\Delta z');
        %title('Average Fisher z-transformed Pearson correlation in ROIs');
        %ylabel('z');

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        %legend(fields(ix), 'interpreter', 'none');
        l = legend({'object update', 'relation update', 'goal update'}, 'interpreter', 'none');
        l.Position = [0.1843 0.7667 0.2037 0.1000];
        xticklabels({'Frontal/Motor (IFG)', 'Dorsal/Parietal', 'Ventral/Temporal'});

         
        orient(gcf, 'landscape');
        print('pdf/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_grouped_GP_components.pdf', '-dpdf');
    
    %
    %
    %% Figure 7: state representations
    %
    %

    case 'plot_gp_CV_state'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=state__nsamples=100_project=0_norm=1_fast=1.mat'));  % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_irrelevant'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=irrelevant__nsamples=100_project=0_norm=1_fast=1.mat'));  % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);


    case 'plot_gp_CV_rois_state_fraction_AAL2_GP_EMPA'
        % plot_gp_CV_rois.m
        
        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_no_project.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 519 725 420]);
        ix = ismember(regressor_names, {'theory', 'state', 'irrelevant', 'PCA'});
        jx = ismember(roi_names, {'LING', 'CAL', 'CUN'});
        h = plot_gp_CV_rois_helper(fs(jx,ix,:), 'signrank', 'median', regressor_names(ix), roi_names(jx), [], [], 10);
        title('Model comparison by ROI');
        ylabel('Fraction significant voxels');

        % Prettifyit
        text(1, 0.75, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([1.5 1.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(2.5, 0.75, 'Early visual', 'fontsize', 12, 'HorizontalAlignment', 'center');
        legend({'EMPA', 'PCA', 'state', 'irrelevant'});
        

    case 'plot_gp_CV_rois_state_fraction_AAL2_GP_EMPA_grouped'
        % plot_gp_CV_rois.m
        
        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped_no_project.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 519 725 420]);
        ix = ismember(regressor_names, {'theory', 'state', 'irrelevant', 'PCA'});
        jx = ismember(roi_names, {'Early visual'});
        h = plot_gp_CV_rois_helper(fs(jx,ix,:), 'signrank', 'median', regressor_names(ix), roi_names(jx), [], [], 10);
        title('Model comparison by ROI group');
        ylabel('Fraction significant voxels');
        xticklabels({'Early visual'});
        legend({'EMPA', 'PCA', 'state', 'irrelevant'});

    otherwise
        assert(false, 'Invalid figure name');
end
