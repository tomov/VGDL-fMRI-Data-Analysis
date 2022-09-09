function show_figure(figure_name)

    figure_scale = 0.7;

    switch figure_name

        case 'example_HRRs'

            rng default;

            D = 100;
            color = gen_embedding(D);
            red = gen_embedding(D);
            blue = gen_embedding(D);
            type = gen_embedding(D);
            avatar = gen_embedding(D);
            missile = gen_embedding(D);


            X(1,:) = bind(type, avatar) + bind(color, blue);
            X(2,:) = bind(type, avatar) + bind(color, red);
            X(3,:) = bind(type, missile) + bind(color, blue);
            X(4,:) = bind(type, missile) + bind(color, red);

            % Normalize
            X = rdivide(X, sqrt(sum(X .^ 2, 2)));

            labels = {'type * avatar + color * blue', ...
                      'type * avatar + color * red', ...
                      'type * missile + color * blue', ...
                      'type * missile + color * red'};
            labels = {'A', 'B', 'C', 'D'};

            figure('pos', [712 152 300*figure_scale 800*figure_scale]);

            subplot(3,1,1);
            imagesc(X);
            title('Example HRRs');
            xlabel('component');
            yticks([1,2,3,4]);
            yticklabels(labels);
            colorbar;

            subplot(3,1,2);
            imagesc(corr(X'));
            xticks([1,2,3,4]);
            yticks([1,2,3,4]);
            xticklabels(labels);
            yticklabels(labels);
            title('Correlation matrix');
            %xtickangle(60);
            %ytickangle(60);
            colorbar;

            subplot(3,1,3);
            imagesc(X * X');
            xticks([1,2,3,4]);
            yticks([1,2,3,4]);
            xticklabels(labels);
            yticklabels(labels);
            title('GP kernel');
            %xtickangle(60);
            %ytickangle(60);
            colorbar;

            print('svg/neuron_revision/figure_neuron_R1_regressors_example_HRRs.svg', '-dsvg');
        
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
            title('t-SNE');
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

            print('svg/neuron_revision/figure_neuron_R1_regressors_theory_HRRs_GLM1_projected.svg', '-dsvg');




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

        case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA__components_and_EMPA__outliers'
            % plot_gp_CV_rois.m
            % see figure3.m

            agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_25M_e1k.mat');
            agg_filename
            load(agg_filename);

            figure('position', [147 521 1045 318]);
            ix = ismember(regressor_names, {'theory', 'sprite', 'interaction', 'termination'});
            cmap = [1 0.7 0.4 0.1]' * [0    0.4470    0.7410] + [0.0 0.3 0.6 0.9]' * [0 0 0];
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
            l = legend({'theory', 'objects', 'relations', 'goals'});
            l.Position = [0.1371 0.7119 0.0842 0.2121];

            orient(gcf, 'landscape');
            print('svg/neuron_revision/plot_gp_CV_rois_fraction_AAL2_GP_EMPA__components_and_EMPA__outliers.svg', '-dsvg');



        case 'plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__components_and_EMPA__outliers'
            % plot_gp_CV_rois.m

            agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.050_atlas=AAL2_GP_EMPA2_grouped_25M_e1k.mat');
            agg_filename
            load(agg_filename);

            figure('position', [147 521 345 318]);
            ix = ismember(regressor_names, {'theory', 'sprite', 'interaction', 'termination'});
            cmap = [1 0.7 0.4 0.1]' * [0    0.4470    0.7410] + [0.0 0.3 0.6 0.9]' * [0 0 0];
            h = plot_gp_CV_rois_helper_boxcharts(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap, 8, [1:3], 1.0, true);
            title('EMPA theory components');
            ylabel('Fraction significant voxels');
            xticklabels({'Frontal/Motor', 'Dorsal/Parietal', 'Ventral/Temporal', 'Early visual'}); 

            % prettify
            yscale = 11.9;
            xscale = 1.25;
            xlim([0.5 * xscale xscale * (length(roi_names) + 0.5)]);
            l = legend({'theory', 'objects', 'relations', 'goals'});
            l.Position = [0.1593 0.7034 0.2481 0.1934];

            print('svg/neuron_revision/plot_gp_CV_rois_fraction_AAL2_GP_EMPA_grouped__components_and_EMPA__outliers.svg', '-dsvg'); 
     

        otherwise
            assert(false, 'Invalid figure name');
    end
end


function c = bind(a, b)
    % Circular convolution
    c = ifft(fft(a) .* fft(b));
end

function embedding = gen_embedding(D)
    embedding = randn(1, D) / sqrt(D);
end
