function show_figure(figure_name)

dqn_pca = true;


switch figure_name

    case 'plot_gp_CV_EMPA'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1.mat')); % this is it !!!!!!!!!!!!!!!!!!!!ontroller five
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    case 'plot_gp_CV_EMPA_table'
        % plot_gp_CV.m

        load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1.mat'));
        filename = bspmview_save_map(vgdl_expt, tmap);
        %ccnl_results_table('AAL2', 'peak', filename, [], [], 0.001, '+/-', 0.05, 20, 3, true, [], 31);
        ccnl_results_table('AAL2', 'vote', filename, [], [], 0.001, '+', 0.05, 20, 3, true, [], 31);

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
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA_25M_e1k.mat');
        if dqn_pca
            agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_neuron.mat';
        else
            agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_25M_e1k.mat');  % !!!! THIS IS IT -- initial submission
        end
        %agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.001_atlas=AAL2_GP_EMPA_25M_e1k.mat');
        agg_filename
        load(agg_filename);

        figure('position', [147 521 1045 218]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});
        %h = plot_gp_CV_rois_helper(fs(:,ix,:), 'ranksum', 'median', regressor_names(ix), roi_names, [], [], 1, [1:2]);
        %h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 2, [1:1]);
        %h = plot_gp_CV_rois_helper(zs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names, [], [], 1, [1:2]);
        h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 8, [1:1], 0.9);
        title('Model comparison by ROI');
        ylabel('Fraction significant voxels');

        % Prettify it
        yscale = 7.2;
        xscale = 1.25;
        text(xscale * 4.5, yscale * 0.075, 'Frontal/Motor', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 7.5 xscale * 7.5], yscale * [0 0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 11.0, yscale * 0.075, 'Dorsal/Parietal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 14.5 xscale * 14.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 16.8, yscale *0.075, 'Ventral/Temporal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 18.5 xscale * 18.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 20.0, yscale * 0.075, 'Early visual', 'fontsize', 10, 'HorizontalAlignment', 'center');
        ylim([0 yscale *0.08]);
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        l = legend({'EMPA', 'DDQN', 'PCA', 'VAE'});
        l.Position = [0.1371 0.6119 0.0842 0.2821];

        orient(gcf, 'landscape');
        if dqn_pca
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_neuron.svg', '-dsvg');
        else
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA.svg', '-dsvg');
        end


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA__outliers'
        % plot_gp_CV_rois.m

        if dqn_pca
            agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_neuron.mat';
        else
            agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_25M_e1k.mat');
        end
        agg_filename
        load(agg_filename);

        figure('position', [147 521 1045 318]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});
        h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 8, [1:1], 0.9, true);
        title('Model comparison by ROI');
        ylabel('Fraction significant voxels');

        % Prettify it
        yscale = 11.9;
        xscale = 1.25;
        text(xscale * 4.5, yscale * 0.075, 'Frontal/Motor', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 7.5 xscale * 7.5], yscale * [0 0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 11.0, yscale * 0.075, 'Dorsal/Parietal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 14.5 xscale * 14.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 16.8, yscale *0.075, 'Ventral/Temporal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 18.5 xscale * 18.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 20.0, yscale * 0.075, 'Early visual', 'fontsize', 10, 'HorizontalAlignment', 'center');
        ylim([0 yscale *0.08]);
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        l = legend({'EMPA', 'DDQN', 'PCA', 'VAE'});
        l.Position = [0.1371 0.7219 0.0842 0.2121];

        orient(gcf, 'landscape');
        if dqn_pca
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA__outliers_neuron.svg', '-dsvg');
        else
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA__outliers.svg', '-dsvg');
        end

        % stats
        g_ROI = []; % ROI
        g_macroROI = []; % ROI group
        g_model = []; % model
        f = fs(:,ix,:);
        ff = [];
        rois = roi_names(:);
        models = regressor_names(ix);
        for i = 1:size(f,1)
            for j = 1:size(f,2)
                ff = [ff; squeeze(f(i,j,:))];
                g_ROI = [g_ROI; repmat(rois(i), [size(f,3) 1])];
                g_model = [g_model; repmat(models(j), [size(f,3) 1])];
            end
        end
        g_macroROI = cell(size(g_ROI));
        g_macroROI(ismember(g_ROI, roi_names(1:7))) = {'frontal'};
        g_macroROI(ismember(g_ROI, roi_names(8:14))) = {'parietal'};
        g_macroROI(ismember(g_ROI, roi_names(15:18))) = {'ventral'};
        g_macroROI(ismember(g_ROI, roi_names(19:21))) = {'visual'};

        [p,tbl,stats] = anovan(ff, {g_macroROI, g_model}, 'model','interaction','varnames',{'macroROI','model'})


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers'
        % plot_gp_CV_rois.m

        if dqn_pca
            agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_neuron.mat';
        else
            agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_25M_e1k.mat');
        end
        agg_filename
        load(agg_filename);

        figure('position', [147 521 345 318]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});
        h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 8, [1:1], 1.0, true);
        title('Model comparison by ROI group');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 

        % prettify
        yscale = 11.9;
        xscale = 1.25;
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        l = legend({'EMPA', 'DDQN', 'PCA', 'VAE'});
        l.Position = [0.1393 0.7334 0.2481 0.1934];

        if dqn_pca
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers_neuron.svg', '-dsvg'); 
        else
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers.svg', '-dsvg'); 
        end
 
        % stats
        % rows = ROIs
        % cols = models
        % reps = subjects
        f = fs(:,ix,:);
        ff = [];
        for i = 1:size(f,1)
            ff = [ff; squeeze(f(i,:,:))'];
        end
        reps = size(fs, 3);
        [p,tbl,stats] = friedman(ff, reps)
        [p,tbl,stats] = anova2(ff, reps)


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA__components__outliers'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_25M_e1k.mat');
        agg_filename
        load(agg_filename);

        figure('position', [147 521 1045 318]);
        ix = ismember(regressor_names, {'sprite', 'interaction', 'termination'});
        cmap = [0.9 0.5 0.2]' * [0    0.4470    0.7410] + [0.1 0.5 0.8]' * [0 0 0];
        h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap, 8, [1:3], 0.9, true);
        title('EMPA theory components');
        ylabel('Fraction significant voxels');

        % Prettify it
        yscale = 11.9;
        xscale = 1.25;
        text(xscale * 4.5, yscale * 0.075, 'Frontal/Motor', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 7.5 xscale * 7.5], yscale * [0 0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 11.0, yscale * 0.075, 'Dorsal/Parietal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 14.5 xscale * 14.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 16.8, yscale *0.075, 'Ventral/Temporal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 18.5 xscale * 18.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 20.0, yscale * 0.075, 'Early visual', 'fontsize', 10, 'HorizontalAlignment', 'center');
        ylim([0 yscale *0.08]);
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        l = legend({'objects', 'relations', 'goals'});
        l.Position = [0.1371 0.7119 0.0842 0.2121];

        orient(gcf, 'landscape');
        print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA__components__outliers.svg', '-dsvg');



    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__components__outliers'
        % plot_gp_CV_rois.m

        agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_25M_e1k.mat');
        agg_filename
        load(agg_filename);

        figure('position', [147 521 345 318]);
        ix = ismember(regressor_names, {'sprite', 'interaction', 'termination'});
        cmap = [0.9 0.5 0.2]' * [0    0.4470    0.7410] + [0.1 0.5 0.8]' * [0 0 0];
        h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap, 8, [1:3], 1.0, true);
        title('EMPA theory components');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 

        % prettify
        yscale = 11.9;
        xscale = 1.25;
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        l = legend({'objects', 'relations', 'goals'});
        l.Position = [0.1593 0.7034 0.2481 0.1934];

        print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__components__outliers.svg', '-dsvg'); 
 


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA__DDQN_layers__outliers'
        % plot_gp_CV_rois.m

        if dqn_pca
            agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_neuron.mat';
        else
            agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_25M_e1k.mat'); % !!!! initial submission
        end
        agg_filename
        load(agg_filename);

        figure('position', [147 521 1045 318]);
        ix = ismember(regressor_names, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
        cmap = [0.9 0.7 0.5 0.3 0.1]' * [0.8500    0.3250    0.0980] + [0.1 0.3 0.5 0.7 0.9]' * [0 0 0];
        h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap, 8, [1:5], 0.9, true);
        title('DDQN layers');
        ylabel('Fraction significant voxels');

        % Prettify it
        yscale = 12.9;
        xscale = 1.25;
        text(xscale * 4.5, yscale * 0.075, 'Frontal/Motor', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 7.5 xscale * 7.5], yscale * [0 0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 11.0, yscale * 0.075, 'Dorsal/Parietal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 14.5 xscale * 14.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 16.8, yscale *0.075, 'Ventral/Temporal', 'fontsize', 10, 'HorizontalAlignment', 'center');
        plot([xscale * 18.5 xscale * 18.5], [0 yscale *0.18], '--', 'color', [0.5 0.5 0.5]);
        text(xscale * 20.0, yscale * 0.075, 'Early visual', 'fontsize', 10, 'HorizontalAlignment', 'center');
        ylim([0 yscale *0.08]);
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        l = legend({'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
        l.Position = [0.1371 0.7119 0.0842 0.2121];

        orient(gcf, 'landscape');
        if dqn_pca
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA__DDQN_layers__outliers__neuron.svg', '-dsvg');
        else
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA__DDQN_layers__outliers.svg', '-dsvg');
        end


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__DDQN_layers__outliers'
        % plot_gp_CV_rois.m

        if dqn_pca
            agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_neuron.mat';
        else
            agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_25M_e1k.mat');
        end
        agg_filename
        load(agg_filename);

        figure('position', [147 521 345 318]);
        ix = ismember(regressor_names, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
        cmap = [0.9 0.7 0.5 0.3 0.1]' * [0.8500    0.3250    0.0980] + [0.1 0.3 0.5 0.7 0.9]' * [0 0 0];
        h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap, 8, [1:5], 1.0, true);
        title('DDQN layers');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 

        % prettify
        yscale = 12.9;
        xscale = 1.25;
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        l = legend({'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
        l.Position = [0.1693 0.7534 0.1881 0.1334];

        if dqn_pca
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__DDQN_layers__outliers_neuron.svg', '-dsvg'); 
        else
            print('svg/figure3/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__DDQN_layers__outliers.svg', '-dsvg'); 
        end
 


    case 'plot_ridge_CV_EMPA'
        % plot_ridge_CV.m

        load(fullfile(get_mat_dir(), 'agg_ridge_CV_us=1_glm=1_model=EMPA_theory_subsample=0_project=1.mat'));
        assert(use_smooth);
        EXPT = vgdl_expt();
        bspmview_wrapper(EXPT, tmap);

    otherwise
        assert(false, 'Invalid figure name');
end
