function show_figure(figure_name)


switch figure_name

    case 'plot_gp_CV_rois_fraction_ungrouped'

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
end
