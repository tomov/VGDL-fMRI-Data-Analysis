close all;
clear all;

EXPT = vgdl_expt();

subjects = 1:1:length(EXPT.subject);

agg_filename = fullfile(get_mat_dir(), sprintf('glm_bic_bms_CV.mat'));
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
nsubjects = length(subjects)


% calc BICs
bic = nan(nROIs, nglmodels, nsubjects); 
lme = nan(nROIs, nglmodels, nsubjects); 

for g = 1:nglmodels
    g
    glmodel = glmodels(g);

    for s = 1:nsubjects
        s
        subj_id = subjects(s);

        %filename = sprintf('glm_predict_glm=%d_subj=%d.mat', glmodel, subj_id);
        filename = sprintf('glm_predict4_glm=%d_subj=%d.mat', glmodel, subj_id);
        load(fullfile(get_mat_dir(2), filename), 'ResMS', 'SPM');

        for m = 1:nROIs
            [N,K] = size(SPM.xX.X);
            K=0; % CV => no free parameters
            % see ccnl_bic.m
            bic(m,g,s) = N*nansum(log(ResMS(roi_masks_flattened{m}))) + K*log(N);
        end
    end
end

lme = -0.5 * bic;

% run BMS
for m = 1:nROIs
    [alpha, exp_r, xp, pxp, bor] = bms(squeeze(lme(m,:,:))'); % col=model, row=subj
    pxp
    pxps(c,:) = pxp;
end

glmodels
table(regions, pxps)



agg_filename

save(agg_filename);
