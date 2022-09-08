function show_figure(figure_name)


switch figure_name

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

        %figure('pos', [64 460 1296*0.5 799*0.5]);
        figure('position', [147 605 1211 134]);

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

            subplot(1, 7, m);
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

            end

            if ismember(m, [7])
                ylim([ax.YLim(1) ax.YLim(2) * 1.3]);
            end
            plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
            plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            if m == 1
                %l = legend(hh, fields, 'interpreter', 'none');
                %l = legend(hh, {'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
                %l.Position = [0.4023 0.1654 0.1289 0.1314];

                %ylabel('z');
                ylabel('\Delta z');
            end
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
        end

        %orient(gcf, 'landscape');
        %print('pdf/plot_PETH_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf', '-bestfit');
        print('svg/figure6/plot_PETH_AAL2_GP_EMPA_GLM_102_GP.svg', '-dsvg');
 

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

        %figure('pos', [64 460 1296*0.5 799*0.5]);
        figure('position', [147 605 1211 134]);

        % only look at the components
        fields = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

        subjs = 1:1:32;

        %cmap = colormap(jet(length(fields)));
        %cmap = colormap(jet(9));
        %cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        cmap = [0.9 0.5 0.2]' * [0    0.4470    0.7410] + [0.1 0.5 0.8]' * [0 0 0];
        t = PETH_dTRs * EXPT.TR; % s

        % loop over masks
        for m = 1:length(mask_filenames)
            disp(mask_name{m});

            subplot(1, 7, m);
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
                            %h = text(t(j) + ix * 0.5, ax.YLim(2) + 0.005 + ix * 0.002, significance(p(j)), 'color', cmap(i,:), 'fontsize', 7, 'HorizontalAlignment', 'center');
                            %set(h,'Rotation',90);
                        end
                    end
                end

            end

            plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);
            plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            if m == 1
                %legend(hh, fields, 'interpreter', 'none');
                %l = legend(hh, {'object update', 'relation update', 'goal update'}, 'interpreter', 'none');
                %l.Position = [0.4162 0.1941 0.1073 0.0588];
                %ylabel('z');
                ylabel('\Delta z');
            end
            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
        end

        %orient(gcf, 'landscape');
        %print('pdf/plot_PETH_components_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf', '-bestfit');
        print('svg/figure6/plot_PETH_components_AAL2_GP_EMPA_GLM_102_GP.svg', '-dsvg');


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
        %figure('pos', [49 329 2143*0.5 610*0.5]);
        figure('position', [147 521 1045 258]);

        ix = 1:nregressors;
        %h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 2, 1, 4);
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, [], cmap, 2.5, 1, 2);
        %title('Average Fisher z-transformed Pearson correlation change');
        ylabel('\Delta z');
        %title('Average Fisher z-transformed Pearson correlation in ROIs');
        %ylabel('z');

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.18, 'Frontal/Motor', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.3], '--', 'color', [0.5 0.5 0.5]);
        text(3.5, 0.18, 'Dorsal/Parietal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([4.5 4.5], [0 0.3], '--', 'color', [0.5 0.5 0.5]);
        text(6, 0.18, 'Ventral/Temporal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        %l = legend(fields(ix), 'interpreter', 'none');
        l = legend({'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
        l.Position = [0.5597 0.7238 0.0807 0.1836];
        l.FontSize = 7;
        ylim([-0.05 0.2]);

        %orient(gcf, 'landscape');
        %print('pdf/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf', '-bestfit');
        print('svg/figure6/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP.svg', '-dsvg');



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
        %cmap = [0.8 0.5 0.2]' * [0    0.4470    0.7410] + [0.2 0.5 0.8]' * [1 1 1];
        cmap = [0.9 0.5 0.2]' * [0    0.4470    0.7410] + [0.1 0.5 0.8]' * [0 0 0];
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
        %figure('pos', [750 464 1457 471]);
        figure('position', [147 521 1045 258]);

        ix = 1:nregressors;
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 5, 1:3, 1.5);
        %title('Average Fisher z-transformed Pearson correlation change');
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
        l.Position = [0.5404 0.6824 0.0988 0.1062];

        %orient(gcf, 'landscape');
        %print('pdf/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP_components.pdf', '-dpdf', '-bestfit');
        print('svg/figure6/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP_components.svg', '-dsvg');



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

    otherwise
        assert(false, 'Invalid figure name');
end
