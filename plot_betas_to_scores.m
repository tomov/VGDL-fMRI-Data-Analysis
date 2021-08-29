

close all;
clear all;

%load('mat/betas_to_scores_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_reg=theory_change_flag.mat');
load('mat/betas_to_scores_glm=102_con=theory_change_flag_Num=1_sphere=10.0mm_reg=theory_change_flag.mat');


%ROI_ix = [1     2     3     5     7    11];
ROI_ix = [1         7    11];
mask_filenames = mask_filenames(ROI_ix);
mask_name = mask_name(ROI_ix);
regions = regions(ROI_ix);
ts = ts(ROI_ix);
betas = betas(ROI_ix);
scores = scores(ROI_ix);

figure('position', [97 451 2190 888]);

clear ps;
clear rhos;

for m = 1:length(mask_filenames)
    subplot(3,5,m);
    hold on;

    [r,p] = corr(scores{m}', ts{m}', 'type', 'spearman');
    ps(m,:) = p;
    rhos(m,:) = r;
    p_corr = 1 - (1 - p)^length(mask_filenames);

    scatter(ts{m}, scores{m});
    lsline;
    text(0.4, 3.0, sprintf('\\rho=%.2f, p_{corr.}=%.2f', r, p_corr), 'fontsize', 15, 'interpreter', 'tex');

    xlabel('t statistic');
    ylabel('expected payout');
    title({regions{m}, mask_name{m}}, 'interpreter', 'none');
end

ps_corr = bonferroni(ps);
table(regions, mask_name', rhos, ps, ps_corr)
