

close all;
clear all;

load('mat/betas_to_scores_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_reg=theory_change_flag.mat');

figure('position', [97 451 2190 888]);

for m = 1:length(mask_filenames)
    subplot(3,5,m);
    hold on;

    [r,p] = corr(wins{m}', betas{m}', 'type', 'spearman');

    scatter(betas{m}, wins{m});
    lsline;
    text(0.4, 3.0, sprintf('\\rho=%.2f, p=%.2f', r, p), 'fontsize', 15, 'interpreter', 'tex');

    xlabel('beta coefficient');
    ylabel('average wins');
    title({regions{m}, mask_name{m}}, 'interpreter', 'none');
end

