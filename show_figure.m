function show_figure(figure_name)


switch figure_name

    case 'plot_gp_CV_rois_fraction_ungrouped'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_ungrouped.mat');
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
        text(14.5, 0.075, 'Early visual', 'fontsize', 12, 'HorizontalAlignment', 'center');
        legend({'EMPA', 'DDQN', 'PCA'});


    case 'plot_gp_CV_rois_fraction_grouped3'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_grouped3.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 522 537 417]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
        h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names);
        title('Model comparison by ROI group');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'});
        legend({'EMPA', 'DDQN', 'PCA'});

    case 'plot_gp_CV_rois_fraction_grouped3_EMPA'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_grouped3.mat');
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


    case 'plot_gp_CV_rois_fraction_grouped3_DDQN'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_grouped3.mat');
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


    case 'plot_PETH_components_ungrouped2_BOLD'
        % plot_PETHs.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_ungrouped2_BOLD.mat')); % !!!!!!!!!!!!
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [      1      2      7     10     11     12     13     14     15]; 

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

    case 'plot_PETHs_bars_components_ungrouped2_BOLD'
        % plot_PETHs_bars.m

        load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_ungrouped2_BOLD.mat')); % !!!!!!!!!!!!
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [      1      2      7     10     11     12     13     14     15]; 

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
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_ungrouped.mat');
        text(1.5, 0.75, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(3.5, 0.75, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([4.5 4.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(6, 0.75, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
        plot([7.5 7.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(8.5, 0.75, 'Early visual', 'fontsize', 12, 'HorizontalAlignment', 'center');
        legend(fields(ix), 'interpreter', 'none');

    otherwise
        assert(false, 'Invalid figure name');
end
