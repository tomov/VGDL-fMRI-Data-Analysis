
clear all;
close all;

%agg_filename= '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/glm_predict_rois_alpha=0.050_neuron.mat';
%agg_filename= '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/glm_predict2_rois_alpha=0.050_neuron.mat';
agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/glm_predict3_rois_alpha=0.050_neuron.mat';

agg_filename

load(agg_filename);


figure('position', [1147 521 1045 418]);
ix = [1 5];
model_names = {'theory', 'sprite', 'interaction','termination','all 3'};
%h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', model_names(ix), regions, [], [], 10);
h = plot_gp_CV_rois_helper(zs(:,ix,:), 'ttest', 'mean', model_names(ix), regions, [], [], 10);
%title('Held out subjects');
title('Leave-one-run-out CV');
ylabel('Fisher z-transformed correlation');
%title('Fraction significant voxels in ROIs');
%ylabel('Fraction significant voxels');
