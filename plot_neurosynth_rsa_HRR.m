% plot permutation test results after neurosynth_rsa_HRR.m


close all;
clear all;

load('mat/neurosynth_rsa_HRR_us=0.mat');

if contains(EXPT.rsadir, '_nosmooth')
    EXPT = vgdl_expt_nosmooth();
else
    EXPT = vgdl_expt();
end



alpha = 0.10; % signifinace level for thresholding (uncorr.) based on permutation tests
idx = find([ROI.p] < alpha);


V = spm_vol('masks/mask.nii');

roi_masks = roi_masks(idx);
ROI = ROI(idx);

rho_map = nan(V.dim);
for r = 1:length(roi_masks)
    roi_mask = roi_masks{r};

    rho_map(roi_mask) = ROI(r).t;
end
bspmview_wrapper(EXPT, rho_map);

