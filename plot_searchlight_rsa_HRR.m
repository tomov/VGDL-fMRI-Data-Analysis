% plot t-test results after neurosynth_rsa_HRR.m with searchlights

close all;
clear all;

load('mat/searchlight_rsa_HRR_us=0_r=6.6667_theory.mat');

if contains(EXPT.rsadir, '_nosmooth')
    EXPT = vgdl_expt_nosmooth();
else
    EXPT = vgdl_expt();
end


alpha = 0.05; % signifinace level for thresholding (uncorr.)
idx = find([ROI.p_s] < alpha);


V = spm_vol('masks/mask.nii');

roi_masks = roi_masks(idx);
ROI = ROI(idx);

map = nan(V.dim);
for r = 1:length(ROI)
    map(ROI(r).cor(1), ROI(r).cor(2), ROI(r).cor(3)) = ROI(r).t_s;
end
bspmview_wrapper(EXPT, map);


