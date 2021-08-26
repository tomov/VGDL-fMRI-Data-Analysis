% plot results of gen_PETHs.m, for paired t-tests

close all;
clear all;

%load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=10.0mm_.mat');
load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=10.0mm__.mat');

figure('pos', [64 421 2282 838]);

% optionally plot theory change flag only
fields(find(strcmp(fields, 'theory_change_flag'))) = [];
%fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
%fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
%fields(find(strcmp(fields, 'termination_change_flag'))) = [];

reg_field = 'theory_change_flag';
nuisance_fields = {'effects', 'avatar_collision_flag', 'block_start', 'block_end', 'instance_start', 'instance_end', 'play_start', 'play_end'};

subjs = 2:2:32;
mask_filenames = mask_filenames(1:2);

cmap = colormap(jet(1 + length(nuisance_fields)));
t = PETH_dTRs * EXPT.TR; % s

% loop over masks
for m = 1:length(mask_filenames)
    disp(mask_name{m});

    for i = 1:length(nuisance_fields)
        subplot(length(mask_filenames), length(nuisance_fields), (m - 1) * length(nuisance_fields) + i);
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
                text(t(j), ax.YLim(2) - 0.02, significance(p(j)), 'color', color, 'fontsize', 17, 'HorizontalAlignment', 'center');
            end
        end

        plot([0 0], ax.YLim, '--', 'color', [0.5 0.5 0.5]);

        if m == 1
            legend(hh, fields, 'interpreter', 'none');
        end

        ylabel('\Delta BOLD');
        xlabel('time (s)');
        title({regions{m}, mask_name{m}, [reg_field, ' vs. ', nuisance_field]}, 'interpreter', 'none');
    end

end

