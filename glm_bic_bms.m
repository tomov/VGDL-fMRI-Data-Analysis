% compare GLMs using BMS
% copied from Exploration

clear all;

EXPT = vgdl_expt();
[~,~,goodRuns,goodSubjects] = vgdl_getSubjectsDirsAndRuns();

% compare keypressGLMs in left M1
%
%{
masks = {'masks/m1_L.nii'};
region = {'Left M1'};
%[masks, region] = get_masks(29, 'DV', true, [], 1);
glms = [5 6];
%}

% compare theory_change_flag GLMs
%
masks = {'masks/ClusterMask_spmT_0001_x=48_y=12_z=32_183voxels.nii', 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii', 'masks/ROI_x=48_y=12_z=32_133voxels_Sphere10.nii', 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii'}';
region = {'theory glm3 clust', 'theory glm3 sphere 6', 'theory glm3 sphere 10', 'theory glm3 sphere 1'}';
glms = [59 60 61 3];

%glms = [62 63 64 21];

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

glms
%table(region, pxps)
table(masks, pxps)


save('glm_bic_bms.mat');

