
clear all;
close all;

agg_filename= '';

agg_filename

load(agg_filename);


figure('position', [1147 521 1045 418]);
h = plot_gp_CV_rois_helper(fs, 'signrank', 'median', glmodels, regions, [], [], 10);
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');
