
close all;
clear all;

EXPT = vgdl_expt();

alpha = 0.05; % significance threshold for individual voxels

%atlas = 'AAL2_grouped';
%atlas = 'AAL2_GP_EMPA_grouped';
atlas = 'AAL2_GP_EMPA2';
%atlas = 'AAL2_GP_EMPA2_grouped';
%atlas = 'AAL2_GP_EMPA';
%atlas = 'AAL2_grouped2';
%atlas = 'Brodmann';
%atlas = 'AAL3v1';
%atlas = 'HarvardOxford-maxprob-thr0';

fasse_ncf = false;
%agg_filename = fullfile(get_mat_dir(fasse_ncf), sprintf('gp_CV_rois_alpha=%.3f_atlas=%s_noNoveltyRule.mat', alpha, atlas));
%agg_filename = fullfile(get_mat_dir(fasse_ncf), sprintf('gp_CV_rois_alpha=%.3f_atlas=%s_vae.mat', alpha, atlas));
%agg_filename = fullfile(get_mat_dir(fasse_ncf), sprintf('gp_CV_rois_alpha=%.3f_atlas=%s_25M_e1k.mat', alpha, atlas));
agg_filename = fullfile(get_mat_dir(fasse_ncf), sprintf('gp_CV_rois_by_game_alpha=%.3f_atlas=%s_neuron.mat', alpha, atlas));
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
    fullfile(get_mat_dir(2), 'fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=EMPA_theory_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=1_parts=123_games=%s.mat'), ...
    fullfile(get_mat_dir(2), 'fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=DQN25M_all_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123_games=%s.mat'), ...
    fullfile(get_mat_dir(2), 'fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=PCA__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123_games=%s.mat'), ...
    fullfile(get_mat_dir(2), 'fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=VAE__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123_games=%s.mat'), ...
};

regressor_names = { ...
    'theory', ...
    'DQN', ...
    'PCA', ...
    'VAE', ...
};



nregressors = length(filename_templates)
assert(length(filename_templates) == nregressors);
subjects = 1:length(EXPT.subject);
nsubjects = length(subjects)
ngames = 6

rs = nan(ngames, nROIs, nregressors, nsubjects); % Pearson r's 
zs = nan(ngames, nROIs, nregressors, nsubjects); % Fisher z-transformed Pearson r's 
fs = nan(ngames, nROIs, nregressors, nsubjects); % fraction significant voxel s
lmes = nan(ngames, nROIs, nregressors, nsubjects); % (overfitted) GP LME
bics = nan(ngames, nROIs, nregressors, nsubjects); % BIC using MLE == LME

for reg = 1:nregressors
    reg
    for s = 1:nsubjects
        s
        subj_id = subjects(s);

        game_names = get_game_names_ordered(subj_id);
        for g = 1:length(game_names)
            game_name = game_names{g};

            filename = sprintf(filename_templates{reg}, subj_id, game_name);
            filename
            load(filename, 'r_CV', 'logmarglik');
            r = mean(r_CV, 1); % Pearson r's for all voxels

            n = EXPT.nTRs * 2; % number of data points per partition (two runs)
            t = r .* sqrt(n - 2) ./ sqrt(1 - r.^2); % T statistic for each voxel https://www.statology.org/p-value-correlation-excel/
            p = 2 * (1 - tcdf(t, n - 2)); % two-sided p value for each voxel 
            %p = 1 - tcdf(t, n - 2); % one-sided p value for each voxel 
            significant = p < alpha; % which voxels are significant

            for m = 1:nROIs
                rs(g,m,reg,s) = mean(r(roi_masks_flattened{m}));
                zs(g,m,reg,s) = mean(atanh(r(roi_masks_flattened{m})));
                fs(g,m,reg,s) = mean(significant(roi_masks_flattened{m}));
                lmes(g,m,reg,s) = sum(logmarglik(roi_masks_flattened{m}));
                bics(g,m,reg,s) = 1 * log(EXPT.nTRs) - 2 * lmes(g,m,reg,s); % k = 1 parameter (sigma), n = # TRs
            end
        end
    end
end

agg_filename

save(agg_filename);
