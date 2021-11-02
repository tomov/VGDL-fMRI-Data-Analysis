close all;
clear all;

EXPT = vgdl_expt();

alpha = 0.01; % significance threshold for individual voxels

[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');

roi_labels = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'}, ...
    {'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Frontal_Mid_2_L', 'Frontal_Mid_2_R'}, ...
    {'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
    {'Precentral_L', 'Precentral_R'}, ...
    {'Postcentral_L', 'Postcentral_R'}, ...
    {'Precuneus_L', 'Precuneus_R'}, ...
    {'Supp_Motor_Area_L', 'Supp_Motor_Area_R'}, ...
    {'Fusiform_L', 'Fusiform_R'}, ...
    {'Occipital_Sup_L', 'Occipital_Sup_R'}, ...
    {'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R'}, ...
    {'Lingual_L', 'Lingual_R'}, ...
    {'Calcarine_L', 'Calcarine_R'}, ...
    {'Cuneus_L', 'Cuneus_R'} ...
};

%roi_labels = { ...
%    {'Frontal_Inf_Tri_L'}, ...
%    {'Frontal_Inf_Oper_L'}, ...
%%    {'Precentral_L'}, ...
%};

roi_names = { ...
    'IFG_Tri', ...
    'IFG_Oper', ...
    'Frontal_Mid', ...
    'Frontal_Sup', ...
    'Precentral', ...
    'Postcentral', ...
    'Precuneus', ...
    'SMA', ...
    'Fusiform', ...
    'Occipital_Sup', ...
    'Occipital_Mid', ...
    'Occipital_Inf', ...
    'Lingual', ...
    'Calcarine', ...
    'Cuneus', ...
};

nROIs = length(roi_labels)
assert(length(roi_names) == nROIs);

for m = 1:nROIs
    roi_masks{m} = ccnl_create_mask(roi_labels{m}, fullfile('masks', [roi_names{m}, '.nii']),  'AAL2', true, 'masks/mask.nii');
    roi_masks_flattened{m} = roi_masks{m}(whole_brain_mask);
end

fasse_ncf = false;

filename_templates = { ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=EMPA_theory_nsamples=100_project=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=EMPA_sprite_nsamples=100_project=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=EMPA_interaction_nsamples=100_project=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=EMPA_termination_nsamples=100_project=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=DQN_all_nsamples=100_project=1_norm=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=DQN_conv1_nsamples=100_project=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=DQN_conv2_nsamples=100_project=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=DQN_conv3_nsamples=100_project=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=DQN_linear1_nsamples=100_project=1_fast=1.mat'), ...
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=DQN_linear2_nsamples=100_project=1_fast=1.mat'), ...
};

regressor_names = { ...
    'theory', ...
    'sprite', ...
    'interaction', ...
    'termination', ...
    'DQN', ...
    'conv1', ...
    'conv2', ...
    'conv3', ...
    'linear1', ...
    'linear2', ...
};

nregressors = length(filename_templates)
assert(length(filename_templates) == nregressors);


subjects = 1:length(EXPT.subject);

nsubjects = length(subjects)


rs = nan(nROIs, nregressors, nsubjects); % Pearson r's 
fs = nan(nROIs, nregressors, nsubjects); % fraction significant voxel s

for reg = 1:nregressors
    reg
    for s = 1:nsubjects
        s
        subj_id = subjects(s);

        filename = sprintf(filename_templates{reg}, subj_id);
        load(filename, 'r_CV');
        r = mean(r_CV, 1); % Pearson r's for all voxels

        n = EXPT.nTRs * 2; % number of data points per partition (two runs)
        t = r .* sqrt(n - 2) ./ sqrt(1 - r.^2); % T statistic for each voxel https://www.statology.org/p-value-correlation-excel/
        p = 2 * (1 - tcdf(t, n - 2)); % p value for each voxel 
        significant = p < alpha; % which boxes are significant

        for m = 1:nROIs
            rs(m,reg,s) = mean(r(roi_masks_flattened{m}));
            fs(m,reg,s) = mean(significant(roi_masks_flattened{m}));
        end
    end
end

filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois.mat');
filename

save(filename);
