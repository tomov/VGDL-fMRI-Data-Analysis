
clear all;
close all;

%agg_filename= '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/glm_predict_rois_alpha=0.050_neuron.mat';
agg_filename= '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/glm_predict2_rois_alpha=0.050_neuron.mat';

agg_filename

load(agg_filename);


figure('position', [1147 521 1045 418]);
%h = plot_gp_CV_rois_helper(fs, 'signrank', 'median', {'theory', 'sprite', 'interaction','termination','all 3'}, regions, [], [], 10);
h = plot_gp_CV_rois_helper(zs, 'ttest', 'mean', {'theory', 'sprite', 'interaction','termination','all 3'}, regions, [], [], 10);
%title('Fraction significant voxels in ROIs');
%ylabel('Fraction significant voxels');
