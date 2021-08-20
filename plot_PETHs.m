% plot results of gen_PETHs.m

close all;
clear all;

load('PETHs.mat');
%load('PETHs_glm21_rois1-7.mat');

figure('pos', [64 421 2282 838]);
cmap = colormap(jet(length(fields)));

% loop over masks
for m = 1:length(mask_filenames)
    disp(mask_name{m});

    subplot(3, 6, m);
    hold on;

    for i = 1:length(fields)
        field = fields{i};
        disp(field)

        t = PETH_dTRs * EXPT.TR; % s
        [sem, me] = wse(activations(m).(field));
        %me = nanmean(activations(m).(field), 1);
        %sem = nanstd(activations(m).(field), 1) / sqrt(size(activations(m).(field), 1)); % TODO wse

        %errorbar(dTRs, m, se);
        hh(i) = plot(t, me, 'color', cmap(i,:));
        h = fill([t flip(t)], [me + sem flip(me - sem)], cmap(i,:));
        set(h, 'facealpha', 0.3, 'edgecolor', 'none');
    end

    if m == 1
        legend(hh, fields, 'interpreter', 'none');
    end
    ylabel('\Delta BOLD');
    xlabel('time (s)');
    title(mask_name{m}, 'interpreter', 'none');
end

