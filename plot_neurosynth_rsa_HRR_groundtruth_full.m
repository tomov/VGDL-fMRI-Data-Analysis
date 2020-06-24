% plot t-test results after neurosynth_rsa_HRR_groundtruth_full.m


close all;
clear all;

%load('mat/neurosynth_rsa_HRR_groundtruth_full_us=0.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_full_us=0_nperms=0.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_full_us=0_nperms=0_dist=euclidean.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_full_us=0_nperms=1000_dist=correlation.mat');
load('mat/neurosynth_rsa_HRR_groundtruth_full_us=0_nperms=1000_dist=correlation.mat');

if contains(EXPT.rsadir, '_nosmooth')
    EXPT = vgdl_expt_nosmooth();
else
    EXPT = vgdl_expt();
end

alpha = 1.3; % signifinace level for thresholding (uncorr.) based on permutation tests
%idx = find([ROI.null_p] < alpha);
idx = find([ROI.null_p] < alpha);

V = spm_vol('masks/mask.nii');

roi_masks = roi_masks(idx);
ROI = ROI(idx);

map = nan(V.dim);
for r = 1:length(roi_masks)
    roi_mask = roi_masks{r};

    %map(roi_mask) = ROI(r).null_z;
    map(roi_mask) = ROI(r).avg_rho;
end
bspmview_wrapper(EXPT, map);






