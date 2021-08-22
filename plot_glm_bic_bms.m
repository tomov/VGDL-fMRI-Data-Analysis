close all;
clear all;

load('mat/glm_bic_bms_controls.mat');

subjs = 1:2:32;

clear pxps;

for m = 1:length(mask_filenames)
    %bic = bics{m}(subjs, [2 3 4 8]);
    %bic = bics{m}(subjs, [1 2 3 4 8]);
    bic = bics{m}(subjs, [1  8]);
    %bic = bics{m}(subjs, :);

    lme = -0.5 * bic;
    [alpha, exp_r, xp, pxp, bor] = bms(lme);
    pxps(m,:) = pxp;

    [sem, me] = wse(bic);

    subplot(3,5,m);
    hold on;
    bar(me);
    errorbar(me, sem, 'MarkerSize', 1);

    xticklabels(glm_names);
end

table(regions, pxps)
