% plot t-test results after neurosynth_rsa_HRR_groundtruth_avg.m


close all;
clear all;

%load('mat/neurosynth_rsa_HRR_groundtruth_avg_us=0.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_avg_us=0_nperms=1000.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_avg_us=0_nperms=1000_dist=euclidean.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_avg_us=0_nperms=1000_dist=correlation_game.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_avg_us=0_nperms=1000_dist=correlation_sprite.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_avg_us=0_nperms=1000_dist=correlation_interaction.mat');
load('mat/neurosynth_rsa_HRR_groundtruth_avg_us=0_nperms=1000_dist=correlation_termination.mat');

if contains(EXPT.rsadir, '_nosmooth')
    EXPT = vgdl_expt_nosmooth();
else
    EXPT = vgdl_expt();
end

alpha = 0.05; % signifinace level for thresholding (uncorr.) based on permutation tests
idx = find([ROI.null_p] < alpha);
%idx = find([ROI.p] < alpha);

V = spm_vol('masks/mask.nii');

roi_masks = roi_masks(idx);
ROI = ROI(idx);

map = nan(V.dim);
for r = 1:length(roi_masks)
    roi_mask = roi_masks{r};

    %map(roi_mask) = ROI(r).null_z;
    map(roi_mask) = ROI(r).avg_rho;
    %map(roi_mask) = ROI(r).t;
end
bspmview_wrapper(EXPT, map);











