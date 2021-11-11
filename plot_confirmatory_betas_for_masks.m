close all;
clear all;

%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=85-51-52-.mat');
%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-.mat');

%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=3-85-51-52-.mat');
%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=6.0mm_cglm=3-85-51-52-.mat');
%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=21-86-82-83-.mat');

%load('mat/confirmatory_betas_for_masks_glm=102_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-.mat');

%load('mat/confirmatory_betas_for_masks_glm=102_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-3-85-51-52-.mat');
%load('mat/confirmatory_betas_for_masks_glm=102_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=103-104-105-.mat');
% subselect ROIs
%ROI_ix = [1     2     3     5     7    11];

%load(fullfile(get_mat_dir(false), 'confirmatory_betas_for_masks_atlas=AAL3v1_cglm=102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-.mat'));
load(fullfile(get_mat_dir(false), 'confirmatory_betas_for_masks_atlas=AAL2_ungrouped2_cglm=102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-102-.mat'));
ROI_ix = 1:length(mask_filenames);

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

    subplot(5,6,m);
    %subplot(1,2,m);
    hold on;
    bar(me);
    errorbar(me, sem, 'o', 'MarkerSize', 1);

    ax = gca;
    for j = 1:size(beta, 2)
        if p(j) <= 0.05
            text(j, ax.YLim(2) - 0.1, significance(p(j)), 'fontsize', 7, 'HorizontalAlignment', 'center'); 
        end
    end

    ylabel('beta coefficient');
    ax.TickLabelInterpreter = 'none';
    xticklabels(confirmatory_regressors);
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',7)
    xtickangle(30);
    xticks(1:length(confirmatory_regressors));
    title({regions{m}, mask_name{m}}, 'interpreter', 'none');

    %p_corr = 1 - (1 - p(1)) ^ length(mask_filenames);
    %text(10, mean([0 ax.YLim(2)]), sprintf('p_{corr.} = %.2e', p_corr), 'fontsize', 17);
    ps(m) = p(1);
end

ps_corr = bonferroni(ps)';
table(regions, mask_name', ps', ps_corr)

regions(ps_corr < 0.05)

ROI_ix = find(ps_corr <= 0.05)
