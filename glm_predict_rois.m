close all;
clear all;

EXPT = vgdl_expt();

test_subjects = 2:2:length(EXPT.subject);

alpha = 0.05; % significance threshold for individual voxels

%agg_filename = fullfile(get_mat_dir(), sprintf('glm_predict_rois_alpha=%.3f_neuron.mat', alpha));
agg_filename = fullfile(get_mat_dir(), sprintf('glm_predict3_rois_alpha=%.3f_neuron.mat', alpha));
agg_filename

%% get masks
%
[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');

load(fullfile(get_mat_dir(false), 'glm_bic_bms_atlas=AAL2_GLM_102_multiplex_with_controls.mat'), 'mask_filenames', 'regions', 'glms');
glmodels = glms;
nROIs = length(mask_filenames);

for m = 1:nROIs
    [roi_mask, ~] = ccnl_load_mask(mask_filenames{m});
    roi_masks_flattened{m} = roi_mask(whole_brain_mask);
end

nglmodels = length(glmodels);
nsubjects = length(test_subjects)


rs = nan(nROIs, nglmodels, nsubjects); % Pearson r's 
zs = nan(nROIs, nglmodels, nsubjects); % Fisher z-transformed Pearson r's 
fs = nan(nROIs, nglmodels, nsubjects); % fraction significant voxel s

for g = 1:nglmodels
    g
    glmodel = glmodels(g);

    for s = 1:nsubjects
        s
        subj_id = test_subjects(s);

        %filename = sprintf('glm_predict_glm=%d_subj=%d.mat', glmodel, subj_id);
        filename = sprintf('glm_predict3_glm=%d_subj=%d.mat', glmodel, subj_id);
        load(fullfile(get_mat_dir(2), filename), 'r');


        % TODO dedupe with gp_CV_rois.m
        r = mean(r, 1); % Pearson r's for all voxels; average across partitions, if any

        n = EXPT.nTRs * 2; % number of data points per partition (two runs)
        t = r .* sqrt(n - 2) ./ sqrt(1 - r.^2); % T statistic for each voxel https://www.statology.org/p-value-correlation-excel/
        p = 2 * (1 - tcdf(t, n - 2)); % two-sided p value for each voxel 
        %p = 1 - tcdf(t, n - 2); % one-sided p value for each voxel 
        significant = p < alpha; % which voxels are significant

        for m = 1:nROIs
            rs(m,g,s) = mean(r(roi_masks_flattened{m}));
            zs(m,g,s) = mean(atanh(r(roi_masks_flattened{m})));
            fs(m,g,s) = mean(significant(roi_masks_flattened{m}));
        end
    end
end

agg_filename

save(agg_filename);
