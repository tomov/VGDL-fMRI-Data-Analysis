function show_figure(figure_name)


switch figure_name

    case 'plot_gp_CV_rois_fraction_ungrouped'

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_ungrouped.mat');
        agg_filename
        load(agg_filename);

        figure('position', [1147 521 1045 418]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
        h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names);
        title('Fraction significant voxels in ROIs');
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
        title('Fraction significant voxels in ROIs');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'});
        legend({'EMPA', 'DDQN', 'PCA'});
end
