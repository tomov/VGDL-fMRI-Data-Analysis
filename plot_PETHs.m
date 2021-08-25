% plot results of gen_PETHs.m

close all;
clear all;

%load('mat/PETHs_glm=21_con=theory_change_flag_Num=3_sphere=4.0mm.mat')
%load('mat/PETHs_tag=tomov2018KL_sphere=10.0mm.mat')
%load('mat/PETHs_tag=hayley2021psi_sphere=10.0mm.mat');
%load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=10.0mm.mat');
load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=6.0mm_.mat');

figure('pos', [64 421 2282 838]);

% optionally plot theory change flag only
fields(find(strcmp(fields, 'theory_change_flag'))) = [];
%fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
%fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
%fields(find(strcmp(fields, 'termination_change_flag'))) = [];

subjs = 2:2:32;

cmap = colormap(jet(length(fields)));

% loop over masks
for m = 1:length(mask_filenames)
    disp(mask_name{m});

    subplot(3, 5, m);
    hold on;

    for i = 1:length(fields)
        field = fields{i};
        disp(field)

        t = PETH_dTRs * EXPT.TR; % s
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
                    text(t(j) + ix * 0.1, ax.YLim(2) - 0.02 - ix * 0.01, significance(p(j)), 'color', cmap(i,:), 'fontsize', 17, 'HorizontalAlignment', 'center');
                end
            end
        end
    end

    if m == 1
        legend(hh, fields, 'interpreter', 'none');
    end
    ylabel('\Delta BOLD');
    xlabel('time (s)');
    title({regions{m}, mask_name{m}}, 'interpreter', 'none');
end

