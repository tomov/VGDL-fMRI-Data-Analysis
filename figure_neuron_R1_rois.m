function show_figure(figure_name)


switch figure_name

    case 'plot_PETH_AAL2_GP_EMPA_GLM_102_GP__R1_rois'
        % plot_PETHs.m

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/PETHs_atlas=Brodmann_what=GP__.mat')
        ROI_ix = [5];

        mask_filenames_ = mask_filenames(ROI_ix);
        mask_name_ = mask_name(ROI_ix);
        regions_ = regions(ROI_ix);
        activations_ = activations(ROI_ix);

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/PETHs_atlas=AAL3v1_neuron_what=GP__.mat')
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [7 8 9];

        mask_filenames_ = [mask_filenames_ mask_filenames(ROI_ix)];
        mask_name_ = [mask_name_ mask_name(ROI_ix)];
        regions_ = [regions_; regions(ROI_ix)];
        activations_ = [activations_ activations(ROI_ix)];

        mask_filenames = mask_filenames_;
        mask_name = mask_name_;
        regions = regions_;
        activations = activations_;

        %figure('pos', [64 460 1296*0.5 799*0.5]);
        figure('position', [147 605 700 134]);

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

            subplot(1, 4, m);
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
        print('svg/neuron_revision/plot_PETH_AAL2_GP_EMPA_GLM_102_GP__R1_rois.svg', '-dsvg');



    case 'plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP__R1_rois'
        % plot_PETHs_bars.m

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/PETHs_atlas=Brodmann_what=GP__.mat')
        ROI_ix = [5];

        mask_filenames_ = mask_filenames(ROI_ix);
        mask_name_ = mask_name(ROI_ix);
        regions_ = regions(ROI_ix);
        activations_ = activations(ROI_ix);

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/PETHs_atlas=AAL3v1_neuron_what=GP__.mat')
        %ROI_ix = 1:length(mask_filenames);
        ROI_ix = [7 8 9];

        mask_filenames_ = [mask_filenames_ mask_filenames(ROI_ix)];
        mask_name_ = [mask_name_ mask_name(ROI_ix)];
        regions_ = [regions_; regions(ROI_ix)];
        activations_ = [activations_ activations(ROI_ix)];

        mask_filenames = mask_filenames_;
        mask_name = mask_name_;
        regions = regions_;
        activations = activations_;

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
        figure('position', [147 521 645 258]);

        ix = 1:nregressors;
        %h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 2, 1, 4);
        h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, [], cmap, 2.5, 1, 2);
        ylabel('\Delta z');

        %l = legend(fields(ix), 'interpreter', 'none');
        l = legend({'theory update', 'interaction', 'avatar interaction', 'new object', 'killed object', 'episode start', 'episode end'}, 'interpreter', 'none');
        l.Position = [0.2597 0.5838 0.0807 0.1836];
        l.FontSize = 7;
        xticklabels({'RSC', 'HC', 'PHC', 'TPO'});
        ylim([-0.03 0.11]);

        %orient(gcf, 'landscape');
        %print('pdf/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP.pdf', '-dpdf', '-bestfit');
        print('svg/neuron_revision/plot_PETH_bars_AAL2_GP_EMPA_GLM_102_GP__R1_rois.svg', '-dsvg');




    otherwise
        assert(false, 'Invalid figure name');
end
