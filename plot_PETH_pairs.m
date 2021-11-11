% plot results of gen_PETHs.m, for paired t-tests

close all;
clear all;

%load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=10.0mm_.mat');
%load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=10.0mm__.mat');

%load('mat/PETHs_glm=102_con=theory_change_flag_odd_Num=1_sphere=10.0mm_.mat');
%load(fullfile(get_mat_dir(false), 'PETHs_glm=102_con=theory_change_flag_odd_Num=1_sphere=10.0mm.mat'));
% subselect ROIs
%ROI_ix = [1     2     3     5     7    11];

load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL3v1_BOLD.mat'));
%ROI_ix = 1:length(mask_filenames);
ROI_ix = [11 12 13 14 18 21 22 23 24];

mask_filenames = mask_filenames(ROI_ix);
mask_name = mask_name(ROI_ix);
regions = regions(ROI_ix);
activations = activations(ROI_ix);

figure('pos', [64 421 2282 838]);

%reg_field = 'theory_change_flag';
nuisance_fields = {'effects', 'avatar_collision_flag', 'new_sprites', 'killed_sprites', 'play_start', 'play_end'};
%reg_field = 'sprite_change_flag';
%reg_field = 'interaction_change_flag';
reg_field = 'termination_change_flag';
%nuisance_fields = {'avatar_collision_flag', 'killed_sprites', 'play_start', 'play_end'};

subjs = 1:1:32;
%mask_filenames = mask_filenames(1:2);

cmap = colormap(jet(1 + length(nuisance_fields)));
t = PETH_dTRs * EXPT.TR; % s

% loop over masks
for m = 1:length(mask_filenames)
    disp(mask_name{m});

    for i = 1:length(nuisance_fields)
        %subplot(length(mask_filenames), length(nuisance_fields), (m - 1) * length(nuisance_fields) + i);
        subplot(length(nuisance_fields), length(mask_filenames), (i - 1) * length(mask_filenames) + m);
        hold on;

        nuisance_field = nuisance_fields{i};

        fields = {reg_field, nuisance_field};

        clear hh;
        for k = 1:length(fields)
            field = fields{k};

            D = activations(m).(field)(subjs,:); % subj x TRs PETH's
            [sem, me] = wse(D);
           
            if strcmp(field, reg_field)
                color = cmap(1,:);
            else
                color = cmap(i + 1,:);
            end
            hh(k) = plot(t, me, 'color', color);
            h = fill([t flip(t)], [me + sem flip(me - sem)], color);
            set(h, 'facealpha', 0.3, 'edgecolor', 'none');
        end

        D1 = activations(m).(reg_field)(subjs,:);
        D2 = activations(m).(nuisance_field)(subjs,:);
        [h,p,ci,stats] = ttest(D1, D2);

        ax = gca;
        xlim([t(1) - 1, t(end) + 1]);
        j_init = find(PETH_dTRs > 0); % ignore baselines
        for j = j_init:length(t)
            if p(j) <= 0.05
                if stats.tstat(j) > 0
                    color = cmap(1,:);
                else
                    color = cmap(i+1,:);
                end
                text(t(j), ax.YLim(2) - 0.02, significance(p(j)), 'color', color, 'fontsize', 7, 'HorizontalAlignment', 'center');
            end
        end

        plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);

        if m == 1
            legend(hh, fields, 'interpreter', 'none');
        end

        ylabel('\Delta BOLD');
        xlabel('time (s)');
        if i == 1
            %title({regions{m}, mask_name{m}, [reg_field, ' vs. ', nuisance_field]}, 'interpreter', 'none');
            title({regions{m}, mask_name{m}, reg_field}, 'interpreter', 'none');
        end
    end

end

