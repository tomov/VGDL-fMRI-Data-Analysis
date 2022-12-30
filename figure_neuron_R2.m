
% In response to neuron reviewer 1, comment about learning

function show_figure(figure_name)

switch figure_name

    case 'corr_theory_HRR_theory_update__subj_1'
        load(fullfile(get_mat_dir(false), 'corr_theory_HRR_theory_update.mat'));

        EXPT = vgdl_expt;
        glmodel = 102;
        subj_id = 1;
        [~, theory_update] = load_GLM_kernel(EXPT, glmodel, subj_id, {'theory_change_flag'}, false, false);
        load(sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', subj_id), 'theory_Xx', 'r_id');
        [coeff,score,latent,tsquared,explained,mu] = pca(theory_Xx);
        
        figure('position', [557 387 386 489]);
        
        subplot(3,1,1)
        hold on;
        plot(theory_Xx(283:2*283,1) * 10);
        plot(theory_update(283:2*283));
        legend({'theory (top PC)', 'theory update'});
        xlabel('time within run (TR)');
        ylabel('signal (a.u.)');
        title('Subject 1');

        subplot(3,1,2)
        bar(cumsum(mean(explaineds(1,1:30), 1)));
        xlabel('Top theory PCs');
        ylabel({'cumulative variance', 'explained (%)'});

        subplot(3,2,5)
        bar(rs(1,1:7));
        hold on;
        xlabel('Top theory PCs');
        ylabel('Pearson''s r');
        %title('                   corr(theory HRR PC, theory update)');

        subplot(3,2,6)
        bar(ps(1,1:7));
        hold on;
        xlabel('Top theory PCs');
        ylabel('p-value');
        %bar(ps(1,:));

        print('svg/neuron_revision/figure_neuron_R2_corr_theory_HRR_theory_update__subj_1.svg', '-dsvg');
        
    
    case 'corr_theory_HRR_theory_update__all_subj'
        load(fullfile(get_mat_dir(false), 'corr_theory_HRR_theory_update.mat'));
        
        zs = atanh(rs);
        zs = zs(:,1:30);
        [sem, me] = wse(zs);
        [h,p,ci,stats] = ttest(zs);
        p

        figure('position', [557 387 386 489]);

        subplot(3,1,1);
        bar(cumsum(mean(explaineds(:,1:30), 1)));
        hold on;
        xlabel('Top theory PCs');
        ylabel({'cumulative variance', 'explained (%)'});
        title('All subjects');

        subplot(3,1,2);
        hold on;
        bar(me);
        h = errorbar(me, sem, '.', 'MarkerSize', 1);
        h.CapSize = 1;
        %ax = gca;
        %for j = 1:size(zs, 2)
        %    if me(j) < 0
        %        y = me(j) - sem(j) - 0.01;
        %    else
        %        y = me(j) + sem(j) + 0.01;
        %    end
        %    h = text(j, y, sprintf('p = %.2f', p(j)), 'fontsize', 6, 'HorizontalAlignment', 'center');
        %    set(h, 'rotation', 60);
        %end
        xlabel('Top theory PCs');
        ylabel('z');
        ylim([-0.08 0.06]);

        subplot(3,1,3);
        bar(p);
        hold on;
        xlabel('Top theory PCs');
        ylabel('p-value');
        %title('                   corr(theory HRR PC, theory update)');

        print('svg/neuron_revision/figure_neuron_R2_corr_theory_HRR_theory_update__all_subj.svg', '-dsvg');


    case 'AAL2_gross'
        [mask_filenames, roi_names, roi_masks] = get_anatomical_masks('AAL2_gross');
        bspmview_wrapper(vgdl_expt, roi_masks{1} + roi_masks{2} * 2 + roi_masks{3} * 3 + roi_masks{4} * 4)


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers__by_game'
        % plot_gp_CV_rois.m

       % agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_25M_e1k.mat');
        agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_by_game_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_project=1_neuron.mat'; % YASSSS
        agg_filename
        load(agg_filename);

        figure('position', [147 521 345 318]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});
        f = squeeze(mean(fs, 1)); % average across games
        h = plot_gp_CV_rois_helper_boxcharts(f(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 8, [1:1], 1.0, true);
        title('Model comparison within games');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 

        % prettify
        yscale = 11.9;
        xscale = 1.25;
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        ylim([0 0.7]);
        l = legend({'EMPA', 'DDQN', 'PCA', 'VAE'});
        l.Position = [0.1393 0.7034 0.2481 0.1934];

        print('svg/neuron_revision/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers__by_game_neuron.svg', '-dsvg'); 

        % stats
        % rows = ROIs
        % cols = models
        % reps = subjects
        f = f(:,ix,:);
        ff = [];
        for i = 1:size(f,1)
            ff = [ff; squeeze(f(i,:,:))'];
        end
        reps = size(fs, 3);
        [p,tbl,stats] = friedman(ff, reps)
        [p,tbl,stats] = anova2(ff, reps)


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers__by_game__stats__DEPRECATED'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_by_game_alpha=0.050_atlas=AAL2_GP_EMPA2_project=1_neuron.mat');
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});

        game_groups = {[1:6], [2,3,4,6], [1,5]};
        game_group_names = {'all', 'more','less'};

        for k = 1:length(game_groups)
            f = squeeze(mean(fs(game_groups{k},:,:,:), 1)); % average across games
            f = f(:,ix,:);

            game_group_names{k}

            % stats
            g_ROI = []; % ROI
            g_macroROI = []; % ROI group
            g_model = []; % model
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
        end
        keyboard

 
    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers__by_game__more_planning'
        % plot_gp_CV_rois.m

        agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_by_game_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_project=1_neuron.mat'; % YASSSS
        agg_filename
        load(agg_filename);

        figure('position', [147 521 345 318]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});
        f = squeeze(mean(fs([2,3,4,6],:,:,:), 1)); % average across games (helper, bait, zelda, lemmings)
        h = plot_gp_CV_rois_helper_boxcharts(f(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 8, [1:1], 1.0, true);
        title('Games with more planning');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 

        % prettify
        yscale = 11.9;
        xscale = 1.25;
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        ylim([0 0.7]);
        l = legend({'EMPA', 'DDQN', 'PCA', 'VAE'});
        l.Position = [0.1393 0.7234 0.2481 0.1934];

        print('svg/neuron_revision/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers__by_game_neuron__more_planning.svg', '-dsvg'); 

        % stats
        % rows = ROIs
        % cols = models
        % reps = subjects
        f = f(:,ix,:);
        ff = [];
        for i = 1:size(f,1)
            ff = [ff; squeeze(f(i,:,:))'];
        end
        reps = size(fs, 3);
        [p,tbl,stats] = friedman(ff, reps)
        [p,tbl,stats] = anova2(ff, reps)

    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers__by_game__less_planning'
        % plot_gp_CV_rois.m

        agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_by_game_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_project=1_neuron.mat'; % YASSSS
        agg_filename
        load(agg_filename);

        figure('position', [147 521 345 318]);
        ix = ismember(regressor_names, {'theory', 'DQN', 'PCA', 'VAE'});
        f = squeeze(mean(fs([1,5],:,:,:), 1)); % average across games (chase,  plaque attack, avoid George)
        h = plot_gp_CV_rois_helper_boxcharts(f(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], [], 8, [1:1], 1.0, true);
        title('Games with less planning');
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 

        % prettify
        yscale = 11.9;
        xscale = 1.25;
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        ylim([0 0.7]);
        l = legend({'EMPA', 'DDQN', 'PCA', 'VAE'});
        l.Position = [0.1193 0.7234 0.2481 0.1934];

        print('svg/neuron_revision/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__outliers__by_game_neuron__less_planning.svg', '-dsvg'); 

        % stats
        % rows = ROIs
        % cols = models
        % reps = subjects
        f = f(:,ix,:);
        ff = [];
        for i = 1:size(f,1)
            ff = [ff; squeeze(f(i,:,:))'];
        end
        reps = size(fs, 3);
        [p,tbl,stats] = friedman(ff, reps)
        [p,tbl,stats] = anova2(ff, reps)


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__more_verses_less_planning'
        % plot_gp_CV_rois.m

        %agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_by_game_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_project=1_neuron.mat'; % YASSSS
        agg_filename = '/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat_from_lab/gp_CV_rois_by_game_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_project=1_neuron.mat'; % YASSSS
        agg_filename
        load(agg_filename);

        figure('position', [147 521 180 318]);
        ix = ismember(regressor_names, {'theory'});
        f1 = squeeze(mean(fs([2,3,4,6],:,ix,:), 1)); % average across games (helper, bait, zelda, lemmings)
        f2 = squeeze(mean(fs([1,5],:,ix,:), 1)); % average across games (chase,  plaque attack, avoid George)
        f1 = reshape(f1, [size(f1, 1), 1, size(f1, 2)]); % add back EMPA dimension
        f2 = reshape(f2, [size(f2, 1), 1, size(f2, 2)]); % add back EMPA dimension
        f = cat(2, f1, f2);
        h = plot_gp_CV_rois_helper_boxcharts(f, 'signrank', 'median', {'more planning', 'less planning'}, roi_names, [], [], 8, [1:1], 1.0, true);
        title({'EMPA by more', 'vs. less planning'});
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 
        xtickangle(50);

        % prettify
        yscale = 11.9;
        xscale = 1.25;
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        ylim([0 0.7]);
        l = legend({'more planning', 'less planning'});
        %l.Position = [0.2093 0.7234 0.2481 0.1934];
        l.Position = [0.5193 0.7834 0.1481 0.0934];

        %set(gcf,'color','none');
        print('svg/neuron_revision/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__more_verses_less_planning.svg', '-dsvg'); 


    case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__approachavoid'
        % plot_gp_CV_rois.m

        agg_filename = '/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat_from_lab/gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_neuron_approachavoid.mat'; %  reviewer 1 approachavoid
        agg_filename
        load(agg_filename);

        figure('position', [147 521 180 318]);
        h = plot_gp_CV_rois_helper_boxcharts(fs, 'signrank', 'median', {'full attributes', 'approach/avoid/neutral'}, roi_names, [], [], 8, [1:1], 1.0, true);
        title({'EMPA by ','sprite encoding'});
        ylabel('Fraction significant voxels');
        xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 
        xtickangle(50);

        % prettify
        yscale = 11.9;
        xscale = 1.25;
        xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
        ylim([0 0.7]);
        l = legend({'full attributes', 'approach/avoid'});
        l.Position = [0.5193 0.7994 0.1481 0.0934];

        %set(gcf,'color','none');
        print('svg/neuron_revision/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__approachavoid.svg', '-dsvg'); 

    otherwise
        assert(false, 'Invalid figure name');
end
