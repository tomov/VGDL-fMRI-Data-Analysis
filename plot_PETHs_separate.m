% plot results of gen_PETHs.m, one in each plot

close all;
clear all;

%load('PETHs.mat');
load('PETHs_glm21_rois1-7.mat');

figure('pos', [64 421 2282 838]);
cmap = colormap(jet(length(fields)));

% loop over masks
for m = 1:length(mask_filenames)
    disp(mask_name{m});

    for i = 1:length(fields)
        subplot(length(mask_filenames), length(fields), (m - 1) * length(fields) + i);
        hold on;

        field = fields{i};
        disp(field)

        t = PETH_dTRs * EXPT.TR; % s
        [sem, me] = wse(activations(m).(field));
        %me = nanmean(activations(m).(field), 1);
        %sem = nanstd(activations(m).(field), 1) / sqrt(size(activations(m).(field), 1)); 

        %errorbar(dTRs, m, se);
        hh = plot(t, me, 'color', cmap(i,:));
        h = fill([t flip(t)], [me + sem flip(me - sem)], cmap(i,:));
        set(h, 'facealpha', 0.3, 'edgecolor', 'none');

        legend(hh, {field}, 'interpreter', 'none');
        ylabel('\Delta BOLD');
        xlabel('time (s)');
        if i == round(length(fields) / 2)
            title(mask_name{m}, 'interpreter', 'none');
        end
    end
end

