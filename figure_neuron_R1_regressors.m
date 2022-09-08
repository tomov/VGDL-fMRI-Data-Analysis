function show_figure(figure_name)

figure_scale = 0.7;

switch figure_name
    
    case 'theory_HRRs'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_Xx');

        figure('pos', [712 152 479*figure_scale 764*figure_scale]);
        imagesc(theory_Xx);
        title('Design matrix');
        xlabel('regressor');
        ylabel('TR');
        colorbar

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_HRRs.svg', '-dsvg');



    case 'theory_kernel'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_kernel');

        figure('pos', [712 152 764*figure_scale 764*figure_scale]);
        imagesc(theory_kernel);
        title('GP kernel');
        xlabel('TR');
        ylabel('TR');
        colorbar

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_kernel.svg', '-dsvg');

    case 'tsne'

        %{
        EXPT = vgdl_expt;
        glmodel = 1;
        subj_id = 1;

        [games, levels] = get_game_for_each_TR(subj_id);

        % Get HRR
        load(sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', subj_id), 'theory_Xx');

        assert(length(games) == size(theory_Xx, 1));
        assert(length(levels) == size(theory_Xx, 1));

        % Run t-SNE
        disp('running tsne...');
        tic
        rng default; % reproducibility
        Y = tsne(theory_Xx);
        toc

        save(fullfile(get_mat_dir(0), 'figure_neuron_R1_regressors_tsne.mat'), 'Y', 'games', 'levels', 'subj_id');
        %}
        

        load(fullfile(get_mat_dir(0), 'figure_neuron_R1_regressors_tsne.mat'), 'Y', 'games', 'levels', 'subj_id');
        proper_games = convert_game_names(games);

        % Plot t-SNE for all TRs
        figure('pos', [712 152 figure_scale*764 figure_scale*764]);
        gscatter(Y(:,1), Y(:,2), proper_games);
        xlabel('dimension 1');
        ylabel('dimension 2');
        title('EMPA theory HRRs: t-SNE');
        legend('Location','southwest');

        print('svg/neuron_revision/figure_neuron_R1_regressors_tsne.svg', '-dsvg');


    case 'tsne_per_game'

        load(fullfile(get_mat_dir(0), 'figure_neuron_R1_regressors_tsne.mat'), 'Y', 'games', 'levels', 'subj_id');

        figure('pos', [712 152 figure_scale*764*3/2+300 figure_scale*764]);

        game_names_ordered = get_game_names_ordered(subj_id);

        proper_game_names_ordered = convert_game_names(game_names_ordered);


        for g = 1:6
            subplot(2, 3, g);
            which_rows = strcmp(games, game_names_ordered{g});

            x = Y(which_rows, 1);
            y = Y(which_rows, 2);
            l = (1:sum(which_rows)) / sum(which_rows) * 9;
            hp = patch([x' NaN], [y' NaN], 0);
            set(hp,'cdata', [l NaN], 'edgecolor','interp','facecolor','none');
            hold on;
            scatter(x,y,15,l,'filled');
            hold off;
            hc = colorbar;
            xlabel('dimension 1');
            ylabel('dimension 2');
            ylabel(hc, 'level', 'FontSize', 11);
            title(proper_game_names_ordered{g});
        end

        print('svg/neuron_revision/figure_neuron_R1_regressors_tsne_per_game.svg', '-dsvg');

    %
    % --------- same but GLM1 projected -------------
    %
    
    case 'theory_HRRs_GLM1_projected'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_Xx');

        EXPT = vgdl_expt;
        glmodel = 1;
        subj = 1;
        [R, K, W] = get_R_K_W(EXPT, glmodel, subj);


        figure('pos', [712 152 479*figure_scale 764*figure_scale]);
        imagesc(R*K*W*theory_Xx);
        title('Design matrix, controlling for game id');
        xlabel('regressor');
        ylabel('TR');

        set(gca,'ColorScale','log');
        h = colorbar;
        %set(h, 'XTickLabel', arrayfun(@(x) sprintf('%.4f', x), exp(h.Ticks), 'UniformOutput', false));

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_HRRs.svg', '-dsvg');




    case 'theory_kernel_GLM1_projected'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_kernel');

        EXPT = vgdl_expt;
        glmodel = 1;
        subj = 1;
        [R, K, W] = get_R_K_W(EXPT, glmodel, subj);
        ker = R*K*W*theory_kernel*W'*K'*R';

        figure('pos', [712 152 764*figure_scale 764*figure_scale]);
        imagesc(ker);
        title('GP kernel, controlling for game id');
        xlabel('TR');
        ylabel('TR');

        set(gca,'ColorScale','log');
        h = colorbar;
        %set(h, 'XTickLabel', arrayfun(@(x) sprintf('%.4f', x), exp(h.Ticks), 'UniformOutput', false));

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_kernel_GLM1_projected.svg', '-dsvg');


    case 'tsne_GLM1_projected'

        %{
        EXPT = vgdl_expt;
        glmodel = 1;
        subj_id = 1;

        [games, levels] = get_game_for_each_TR(subj_id);

        % Get HRR
        load(sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', subj_id), 'theory_Xx');

        assert(length(games) == size(theory_Xx, 1));
        assert(length(levels) == size(theory_Xx, 1));

        % Run t-SNE
        disp('running tsne...');
        tic
        rng default; % reproducibility
        [R, K, W] = get_R_K_W(EXPT, glmodel, subj_id);
        Y = tsne(R * K * W * theory_Xx);
        toc

        save(fullfile(get_mat_dir(0), 'figure_neuron_R1_regressors_tsne_GLM1_projected.mat'), 'Y', 'games', 'levels', 'subj_id');
        %}
        

        load(fullfile(get_mat_dir(0), 'figure_neuron_R1_regressors_tsne_GLM1_projected.mat'), 'Y', 'games', 'levels', 'subj_id');
        proper_games = convert_game_names(games);

        % Plot t-SNE for all TRs
        figure('pos', [712 152 figure_scale*764 figure_scale*764]);
        gscatter(Y(:,1), Y(:,2), proper_games);
        xlabel('dimension 1');
        ylabel('dimension 2');
        title('t-SNE, controlling for game id');
        legend('Location','northeast');

        print('svg/neuron_revision/figure_neuron_R1_regressors_tsne_GLM1_projected.svg', '-dsvg');


    case 'tsne_per_game_GLM1_projected'

        load(fullfile(get_mat_dir(0), 'figure_neuron_R1_regressors_tsne_GLM1_projected.mat'), 'Y', 'games', 'levels', 'subj_id');

        figure('pos', [712 152 figure_scale*764*3/2+300 figure_scale*764]);

        game_names_ordered = get_game_names_ordered(subj_id);

        proper_game_names_ordered = convert_game_names(game_names_ordered);


        for g = 1:6
            subplot(2, 3, g);
            which_rows = strcmp(games, game_names_ordered{g});

            x = Y(which_rows, 1);
            y = Y(which_rows, 2);
            l = (1:sum(which_rows)) / sum(which_rows) * 9;
            hp = patch([x' NaN], [y' NaN], 0);
            set(hp,'cdata', [l NaN], 'edgecolor','interp','facecolor','none');
            hold on;
            scatter(x,y,15,l,'filled');
            hold off;
            hc = colorbar;
            xlabel('dimension 1');
            ylabel('dimension 2');
            ylabel(hc, 'level', 'FontSize', 11);
            title(proper_game_names_ordered{g});
        end

        print('svg/neuron_revision/figure_neuron_R1_regressors_tsne_per_game_GLM1_projected.svg', '-dsvg');


    otherwise
        assert(false, 'Invalid figure name');
end
