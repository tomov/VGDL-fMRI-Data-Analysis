close all;
clear all;

%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=85-51-52-.mat');
load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-21-.mat');

%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=3-85-51-52-.mat');
%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=6.0mm_cglm=3-85-51-52-.mat');
%load('mat/confirmatory_betas_for_masks_glm=21_con=theory_change_flag_Num=1_sphere=10.0mm_cglm=21-86-82-83-.mat');

subjs = 2:2:32;

figure('position', [97 451 2190 888]);

for m = 1:length(mask_filenames)
%for m = 1:2
    beta = betas{m}(subjs, :)
    [sem, me] = wse(beta);
    [h,p,ci,stats] = ttest(beta);

    subplot(3,5,m);
    %subplot(1,2,m);
    hold on;
    bar(me);
    errorbar(me, sem, 'o', 'MarkerSize', 1);

    ax = gca;
    for j = 1:size(beta, 2)
        if p(j) <= 0.05
            text(j, ax.YLim(2) - 0.1, significance(p(j)), 'fontsize', 17, 'HorizontalAlignment', 'center'); 
        end
    end

    ylabel('beta coefficient');
    ax.TickLabelInterpreter = 'none';
    xticklabels(confirmatory_regressors);
    xtickangle(30);
    xticks(1:length(confirmatory_regressors));
    title({regions{m}, mask_name{m}}, 'interpreter', 'none');
end

