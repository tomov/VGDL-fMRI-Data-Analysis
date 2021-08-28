close all;
clear all;

%load('mat/glm_bic_bms_controls.mat');

%load('mat/glm_bic_bms_nocontrols.mat');
%glm_ix = [1 2 3 4 8];

%load('mat/glm_bic_bms_single_controls.mat');
%load('mat/glm_bic_bms_glm=102_con=theory_change_flag_Num=1_sphere=10.0mm_multiplex.mat');
load('mat/glm_bic_bms_glm=102_con=theory_change_flag_Num=1_sphere=10.0mm_multiplex_w_controls.mat');
glm_ix = [1 2 3 4 5];
%ROI_ix = [1     2     3     5     7    11];
%mask_filenames = mask_filenames(ROI_ix);
%regions = regions(ROI_ix);
%bics = bics(ROI_ix);

%subjs = 2:2:32;
subjs = 2:2:32;
tick_labels = {};
for j = 1:length(glm_ix)
    tick_labels{j} = sprintf('GLM %d: %s', glms(glm_ix(j)), glm_names{glm_ix(j)});
end

clear pxps;

figure('position', [97 451 2190 888]);


for m = 1:length(mask_filenames)
    %bic = bics{m}(subjs, [2 3 4 8]);
    %bic = bics{m}(subjs, [1 2 3 4 8]);
    [~, mask_name{m}, ~] = fileparts(mask_filenames{m});
    bic = bics{m}(subjs, glm_ix);
    %bic = bics{m}(subjs, :);

    lme = -0.5 * bic;
    [alpha, exp_r, xp, pxp, bor] = bms(lme);
    pxps(m,:) = pxp;

    [sem, me] = wse(bic);

    subplot(3,5,m);
    hold on;
    bar(me - me(1));
    errorbar(me - me(1), sem, 'o', 'MarkerSize', 1);

    xlabel('GLM');
    ylabel('\Delta BIC');
    %xticklabels(glm_names(glm_ix));
    xticklabels(tick_labels);
    xtickangle(30);
    xticks(1:length(glm_ix));
    title({regions{m}, mask_name{m}}, 'interpreter', 'none');
end

subjs
glm_names(glm_ix)
table(regions, pxps)
