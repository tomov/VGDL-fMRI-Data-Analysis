% compare GLMs using BMS
% copied from Exploration

clear all;

EXPT = vgdl_expt();
[subjects,~,goodRuns,goodSubjects] = vgdl_getSubjectsDirsAndRuns();

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
%{
masks = {'masks/ClusterMask_spmT_0001_x=48_y=12_z=32_183voxels.nii', 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii', 'masks/ROI_x=48_y=12_z=32_133voxels_Sphere10.nii', 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii'}';
region = {'theory glm3 clust', 'theory glm3 sphere 6', 'theory glm3 sphere 10', 'theory glm3 sphere 1'}';
glms = [59 60 61 3];
%}
%glms = [62 63 64 21];

glmodel = 102;
contrast = 'theory_change_flag';
Num = 1;
sphere = 10;
[mask_filenames, regions] = get_masks_from_contrast(glmodel, contrast, true, [], Num, sphere);
mask_filenames = mask_filenames';

% multiplexing
%glms = [21 86 82 83 97 98 88 68];  % with control regressors
%glms = [3 85 51 52 95 96 87 53]; % without controls
%glm_names = {'th', 's', 'i', 't', 's+i', 's+t', 'i+t', 's+i+t'};

% single signal
%glms = [86 82 83 93 94 75 21]; %with control regressors
%glms = [85 51 52 91 92 67 3]; % without control
%glm_names = {'s', 'i', 't', 's|i', 's|t', 'i|t', 'th'};
%filename = 'mat/glm_bic_bms_single_controls.mat';

% multiplexing, again
%glms = [102 103 104 105 106]; % with controls
glms = [3 85 51 52 53]; % w/o controls
glm_names = {'th', 's', 'i', 't', 's+i+t'};
filename = sprintf('mat/glm_bic_bms_glm=%d_con=%s_Num=%d_sphere=%.1fmm_multiplex.mat', glmodel, contrast, Num, sphere);;

for c = 1:length(mask_filenames)
    mask_filename = mask_filenames{c};

    mask_filename

    lmes = [];
    bics{c} = [];
    for i = 1:length(glms)
        glmodel = glms(i);
        bic = ccnl_bic(EXPT, glmodel, mask_filename, subjects);
        bics{c} = [bics{c} bic];
    end

    lme = -0.5 * bics{c};
    [alpha, exp_r, xp, pxp, bor] = bms(lme);

    pxp

    pxps(c,:) = pxp;
end

glms
table(regions, pxps)
%table(mask_filenames, pxps)


save(filename);

