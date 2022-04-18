function figure4(figure_name)


switch figure_name 


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

        %figure('position', [97 451 2190 888]);
        figure('position', [409 4 813 831]);

        %confirmatory_regressors = confirmatory_regressors(end-2:end);
        ps = [];

        for m = 1:length(mask_filenames)
        %for m = 1:2
            beta = betas{m}(subjs, :);
            %beta = beta(:,end-2:end);
            beta

            [sem, me] = wse(beta);
            [h,p,ci,stats] = ttest(beta);

            subplot(6,3,m);
            %subplot(1,2,m);
            hold on;
            bar(me);
            h = errorbar(me, sem, '.', 'MarkerSize', 1);
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
            title(regions{m}, 'interpreter', 'none', 'fontsize', 7);
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
        print('svg/figure4/plot_confirmatory_betas_for_masks_AAL2_GLM_102.svg', '-dsvg');


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

        %figure('pos', [64 348 2110*0.5 911*0.5]);
        figure('position', [147 605 1411 134]);

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

            %subplot(3, 4, m + (floor((m - 1) / 3)));
            subplot(1, 9, m);
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

            end 

            xlabel('time (s)');
            title(regions{m}, 'interpreter', 'none');
            %ylim([-0.2 0.5]); 

            if m <= 6
                ylim([ax.YLim(1) ax.YLim(2) * 1.3]);
            end
            plot(ax.XLim, [0 0], '--', 'color', [0.5 0.5 0.5]);
            plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);

            if m == length(mask_filenames)
                %l = legend(hh, {'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
                %l.Position = [0.9075 0.0703 0.0965 0.8788];
            end
            if m == 1
                if exist('what', 'var') && strcmp(what, 'GP')
                    ylabel('\Delta z');
                else
                    ylabel('\Delta BOLD');
                end 
            end
        end

        %orient(gcf, 'landscape');
        %print('pdf/plot_PETH_AAL2_GLM_102_BOLD.pdf', '-dpdf', '-fillpage');
        print('svg/figure4/plot_PETH_AAL2_GLM_102_BOLD.svg', '-dsvg');


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
        %figure('pos', [49 329 2143 610]);
        figure('position', [147 521 1045 218]);

        ix = 1:nregressors;
        %h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 12, [1:1], 0.4);
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, [], cmap, 10, [1:1], 0.5);
        if exist('what', 'var') && strcmp(what, 'GP')
            title('Average Fisher z-transformed Pearson correlation change');
            ylabel('\Delta z');
        else
            title('Average BOLD change in ROIs');
            ylabel('\Delta BOLD');
        end

        % Prettyfy it 
        % specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
        text(1.5, 0.55, 'Frontal/Motor', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([2.5 2.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(3.00, 0.55, 'Dorsal/Parietal', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([3.5 3.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(5.5, 0.55, 'Ventral/Temporal', 'fontsize', 9, 'HorizontalAlignment', 'center');
        plot([7.5 7.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
        text(8.5, 0.55, 'Early visual', 'fontsize', 9, 'HorizontalAlignment', 'center');
        %legend(fields(ix), 'interpreter', 'none');
        l = legend({'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
        l.Position = [0.6830 0.6410 0.0742 0.2344];
        l.FontSize = 6;
        ylim([-0.1 0.6]);

        % fix the significant stars
        %H = findobj(gcf);
        %H(15).Position=[9.3429 0.1004 0];

        %orient(gcf, 'landscape');
        print('svg/figure4/plot_PETHs_bars_AAL2_GLM_102_BOLD.svg', '-dsvg');



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


    otherwise
        assert(false, 'Invalid figure name');
end
