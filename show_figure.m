function show_figure(figure_name)


switch figure_name

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

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_all_nsamples=100_project=1_norm=1_fast=1.mat')); % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_PCA'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=PCA_all_nsamples=100_project=1_norm=1_fast=1.mat')); % !!!
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 521 1045 418]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
        h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names);
        title('Model comparison by ROI');
        ylabel('Fraction significant voxels');

        % Prettify it
        text(3.5, 0.075, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([6.5 6.5], [0 0.08], '--', 'color', [0.5 0.5 0.5]);
        text(8.5, 0.075, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([10.5 10.5], [0 0.08], '--', 'color', [0.5 0.5 0.5]);
        text(12, 0.075, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([13.5 13.5], [0 0.08], '--', 'color', [0.5 0.5 0.5]);
        text(15, 0.075, 'Early visual', 'fontsize', 12, 'HorizontalAlignment', 'center');
        legend({'EMPA', 'DDQN', 'PCA'});


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 522 537 417]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
        h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names);
        title('Model comparison by ROI group');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'});
        legend({'EMPA', 'DDQN', 'PCA'});

    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped_EMPA_components'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 522 537 417]);
        ix = ismember(regressor_names, {'theory', 'sprite', 'interaction', 'termination'});
        cmap = [1 0.8 0.6 0.4]' * [0    0.4470    0.7410];
        plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap); %colormap(winter(3)));
        title('EMPA components');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'});
        legend({'theory', 'objects', 'relations', 'goals'});


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped_DDQN_layers'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 522 537 417]);
        ix = ismember(regressor_names, {'DQN', 'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
        cmap = [1 0.9 0.8 0.7 0.6 0.5]' * [0.8500    0.3250    0.0980];
        plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap); %colormap(autumn(5)));
        title('DDQN layers');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'});
        legend({'all layers', 'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});


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


    case 'plot_confirmatory_betas_for_masks_AAL2_GLM_102_components'
        % plot_confirmatory_betas_for_masks.m

        load(fullfile(get_mat_dir(false), 'confirmatory_betas_for_masks_atlas=AAL2_GLM_102_cglm=103-104-105-.mat')); % !!!!!!!!!
        ROI_ix = [      1      2      8         11     12     13     14     15  16]; 

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


    case 'plot_PETH_components_AAL2_GLM_102_BOLD'
        % plot_PETHs.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GLM_102_BOLD.mat')); % !!!!!!!!!!!!
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [      1      2      8     11     12     13     14     15   16]; 

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        figure('pos', [64 348 1564 911]);

        % optionally plot theory change flag only
        fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        %fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        %fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        %fields(find(strcmp(fields, 'termination_change_flag'))) = [];
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
                            text(t(j) + ix * 0.1, ax.YLim(2) - 0.02 - ix * 0.02, significance(p(j)), 'color', cmap(i,:), 'fontsize', 7, 'HorizontalAlignment', 'center');
                        end
                    end
                end

                plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
                plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            end

            if m == length(mask_filenames)
                legend(hh, fields, 'interpreter', 'none');
            end
            if exist('what', 'var') && strcmp(what, 'GP')
                ylabel('\Delta z');
            else
                ylabel('\Delta BOLD');
            end
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
        end



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
        fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        %fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        %fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        %fields(find(strcmp(fields, 'termination_change_flag'))) = [];
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
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 5, 1:3);
        if exist('what', 'var') && strcmp(what, 'GP')
            title('Average Fisher z-transformed Pearson correlation change in ROIs');
            ylabel('\Delta z');
        else
            title('Average BOLD change in ROIs');
            ylabel('\Delta BOLD');
        end

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.75, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(3, 0.75, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([3.5 3.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(5.5, 0.75, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([7.5 7.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(8.5, 0.75, 'Early visual', 'fontsize', 12, 'HorizontalAlignment', 'center');
        legend(fields(ix), 'interpreter', 'none');


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
        cmap = [ ...
         0         0    1.0000; ...
         0    0.3333    1.0000; ...
         0    0.6667    1.0000];
        cmap = [cmap; 0.7294    0.3333    0.8275];
        glm_names = {'theory updates', 'object updates', 'relation updates', 'goal updates', 'object, relation, goal updates'};
        ix = [2,3,4,5];
        h = plot_gp_CV_rois_helper(bs(:,ix,:), 'ttest', 'mean', glm_names(ix), regions, 0, cmap, 500000, []);
        ylabel('\Delta BIC relative to theory updates');
        title('GLM model comparison');

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.75 * 50000, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.8 * 50000], '--', 'color', [0.5 0.5 0.5]);
        text(3, 0.75 * 50000, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([3.5 3.5], [0 0.8 * 50000], '--', 'color', [0.5 0.5 0.5]);
        text(5, 0.75 * 50000, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        legend(glm_names(ix),  'interpreter', 'none');


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

        figure('pos', [64 460 1296 799]);

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
                            text(t(j) + ix * 0.1, ax.YLim(2) - 0.02 - ix * 0.01, significance(p(j)), 'color', cmap(i,:), 'fontsize', 7, 'HorizontalAlignment', 'center');
                        end
                    end
                end

                plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
                plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            end

            if m == 1
                legend(hh, fields, 'interpreter', 'none');
            end
            %ylabel('z');
            ylabel('\Delta z');
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
        end

 

    case 'plot_PETH_components_AAL2_GP_EMPA_GLM_102_GP'
        % plot_PETHs.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP.mat')); % !!!!!!!!!!!!
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_CV_no_baseline.mat')); 
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        figure('pos', [64 460 1296 799]);

        % optionally plot theory change flag only
        fields(find(strcmp(fields, 'theory_change_flag'))) = [];
        %fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
        %fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
        %fields(find(strcmp(fields, 'termination_change_flag'))) = [];
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
                            text(t(j) + ix * 0.1, ax.YLim(2) - 0.02 - ix * 0.01, significance(p(j)), 'color', cmap(i,:), 'fontsize', 7, 'HorizontalAlignment', 'center');
                        end
                    end
                end

                plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
                plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            end

            if m == 1
                legend(hh, fields, 'interpreter', 'none');
            end
            ylabel('\Delta z');
            %ylabel('z');
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
        end


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
        figure('pos', [49 329 2143 610]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 5, 1);
        title('Average Fisher z-transformed Pearson correlation change in ROIs');
        ylabel('\Delta z');
        %title('Average Fisher z-transformed Pearson correlation in ROIs');
        %ylabel('z');

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.25, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.3], '--', 'color', [0.5 0.5 0.5]);
        text(3.5, 0.25, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([4.5 4.5], [0 0.3], '--', 'color', [0.5 0.5 0.5]);
        text(6, 0.25, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        legend(fields(ix), 'interpreter', 'none');


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
        figure('pos', [49 272 936 667]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 5, 1);
        title('Average Fisher z-transformed Pearson correlation change in ROIs');
        ylabel('\Delta z');
        %title('Average Fisher z-transformed Pearson correlation in ROIs');
        %ylabel('z');

        % Prettyfy it 
        legend(fields(ix), 'interpreter', 'none');
        xticklabels({'Frontal/Motor (IFG)', 'Dorsal/Parietal', 'Ventral/Temporal'});
         

         

    case 'plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP_components'
        % plot_PETHs_bars.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP.mat')); % !!!!!!!!!!!!
        %load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_CV_no_baseline.mat')); 
        ROI_ix = 1:length(mask_filenames);

        mask_filenames = mask_filenames(ROI_ix);
        mask_name = mask_name(ROI_ix);
        regions = regions(ROI_ix);
        activations = activations(ROI_ix);

        % optionally plot theory change flag only
        fields = {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

        subjs = 1:1:32;

        nROIs = length(mask_filenames);
        nregressors = length(fields);
        nsubjects = length(subjs);

        as = nan(nROIs,nregressors,nsubjects);

        cmap = [1 0.8 0.6 0.4]' * [0    0.4470    0.7410];
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
        figure('pos', [49 329 2143 610]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 5);
        title('Average Fisher z-transformed Pearson correlation change in ROIs');
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
        legend(fields(ix), 'interpreter', 'none');

    
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
