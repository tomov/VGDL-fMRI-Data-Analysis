close all;
clear all;

EXPT = vgdl_expt();

alpha = 0.01; % significance threshold for individual voxels

atlas = 'AAL2_grouped';
%atlas = 'AAL3v1';
%atlas = 'HarvardOxford-maxprob-thr0';

fasse_ncf = false;
agg_filename = fullfile(get_mat_dir(fasse_ncf), sprintf('gp_CV_rois_alpha=%.3f_atlas=%s.mat', alpha, atlas));
agg_filename

%% get masks
%
[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');

[mask_filenames, roi_names, roi_masks] = get_anatomical_masks(atlas);
roi_names = roi_names';
nROIs = length(mask_filenames)

for m = 1:nROIs
    roi_masks_flattened{m} = roi_masks{m}(whole_brain_mask);
end


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
    fullfile(get_mat_dir(fasse_ncf), 'fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=PCA_all_nsamples=100_project=1_norm=1_fast=1.mat'), ...
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
    'PCA', ...
};

nregressors = length(filename_templates)
assert(length(filename_templates) == nregressors);


subjects = 1:length(EXPT.subject);

nsubjects = length(subjects)


rs = nan(nROIs, nregressors, nsubjects); % Pearson r's 
zs = nan(nROIs, nregressors, nsubjects); % Fisher z-transformed Pearson r's 
fs = nan(nROIs, nregressors, nsubjects); % fraction significant voxel s
lmes = nan(nROIs, nregressors, nsubjects); % (overfitted) GP LME
bics = nan(nROIs, nregressors, nsubjects); % BIC using MLE == LME

for reg = 1:nregressors
    reg
    for s = 1:nsubjects
        s
        subj_id = subjects(s);

        filename = sprintf(filename_templates{reg}, subj_id);
        filename
        load(filename, 'r_CV', 'logmarglik');
        r = mean(r_CV, 1); % Pearson r's for all voxels

        n = EXPT.nTRs * 2; % number of data points per partition (two runs)
        t = r .* sqrt(n - 2) ./ sqrt(1 - r.^2); % T statistic for each voxel https://www.statology.org/p-value-correlation-excel/
        p = 2 * (1 - tcdf(t, n - 2)); % p value for each voxel 
        significant = p < alpha; % which boxes are significant

        for m = 1:nROIs
            rs(m,reg,s) = mean(r(roi_masks_flattened{m}));
            zs(m,reg,s) = mean(atanh(r(roi_masks_flattened{m})));
            fs(m,reg,s) = mean(significant(roi_masks_flattened{m}));
            lmes(m,reg,s) = sum(logmarglik(roi_masks_flattened{m}));
            bics(m,reg,s) = 1 * log(EXPT.nTRs) - 2 * lmes(m,reg,s); % k = 1 parameter (sigma), n = # TRs
        end
    end
end

agg_filename

save(agg_filename);
