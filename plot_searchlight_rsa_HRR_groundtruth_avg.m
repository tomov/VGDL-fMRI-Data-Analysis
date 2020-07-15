% plot t-test results after neurosynth_rsa_HRR_groundtruth_avg.m for searchlights


close all;
clear all;

load('mat/searchlight_rsa_HRR_groundtruth_avg_us=0_r=4.0000_nperms=0_dist=correlation_game.mat');

if contains(EXPT.rsadir, '_nosmooth')
    EXPT = vgdl_expt_nosmooth();
else
    EXPT = vgdl_expt();
end

alpha = 0.05; % signifinace level for thresholding (uncorr.) based on permutation tests
idx = find([ROI.p] < alpha);

V = spm_vol('masks/mask.nii');

roi_masks = roi_masks(idx);
ROI = ROI(idx);

map = nan(V.dim);
for r = 1:length(ROI)
    %map(roi_mask) = ROI(r).avg_rho;
    map(ROI(r).cor(1), ROI(r).cor(2), ROI(r).cor(3)) = ROI(r).avg_rho;
    %map(roi_mask) = ROI(r).t;
end
bspmview_wrapper(EXPT, map);











