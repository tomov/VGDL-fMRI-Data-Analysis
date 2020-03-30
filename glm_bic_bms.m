% compare GLMs using BMS
% copied from Exploration

clear all;

EXPT = vgdl_expt();
[~,~,goodRuns,goodSubjects] = vgdl_getSubjectsDirsAndRuns();

% compare GLMs in left M1
%
masks = {'masks/left_M1.nii'};
region = {'Left M1'};
%[masks, region] = get_masks(29, 'DV', true, [], 1);
glms = [5 6];

for c = 1:length(masks)
    mask = masks{c};

    mask

    lmes = [];
    bics{c} = [];
    for i = 1:length(glms)
        glmodel = glms(i);
        bic = ccnl_bic(EXPT, glmodel, mask, goodSubjects);
        bics{c} = [bics{c} bic];
    end

    lme = -0.5 * bics{c};
    [alpha, exp_r, xp, pxp, bor] = bms(lme);

    pxp

    pxps(c,:) = pxp;
end

table(region, pxps)


save('glm_bic_bms.mat');

