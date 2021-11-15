clear all;
close all;


fasse_ncf = false;
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=HarvardOxford-maxprob-thr0.mat');
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL3v1.mat');

%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2.mat'); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_grouped.mat');
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_grouped2.mat');
agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_grouped.mat');
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=Brodmann.mat');

%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010.mat');
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.001.mat');
agg_filename

load(agg_filename);

%
%% BICs
%

%{
figure('position', [673 90 1519 849]);
ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
h = plot_gp_CV_rois_helper(bics(:,ix,:) - bics(:,end,:), 'ttest', 'mean', regressor_names(ix), roi_names);
title('BICs in ROIs');
ylabel('\Delta BIC');
%}


%
%% fraction significant voxels
%

figure('position', [1147 521 1045 418]);
ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names);
%h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, alpha);
%plot_gp_CV_rois_helper(fs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names);
%ylim([0 0.1]);
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');

% Prettify it
% specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
%{
text(3.5, 0.075, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
plot([6.5 6.5], [0 0.08], '--', 'color', [0.5 0.5 0.5]);
text(8.5, 0.075, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
plot([10.5 10.5], [0 0.08], '--', 'color', [0.5 0.5 0.5]);
text(12, 0.075, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
plot([13.5 13.5], [0 0.08], '--', 'color', [0.5 0.5 0.5]);
text(14.5, 0.075, 'Early visual', 'fontsize', 12, 'HorizontalAlignment', 'center');
%}
legend({'EMPA', 'DDQN', 'PCA'});




figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'theory', 'sprite', 'interaction', 'termination'});
cmap = [1 0.8 0.6 0.4]' * h(1).FaceColor;
plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap); %colormap(winter(3)));
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');

figure('position', [173 390 1519 849]);
ix = ismember(regressor_names, {'DQN', 'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
cmap = [1 0.9 0.8 0.7 0.6 0.5]' * h(2).FaceColor;
plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, [], cmap); %colormap(autumn(5)));
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');


%{
figure;
plot_gp_CV_rois_helper(fs, 'signrank', regressor_names, roi_names);
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');
%}

%% Pearson correlations
%

%{
figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
h = plot_gp_CV_rois_helper(zs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names, 0);
ylim([0 0.04])
title('Fisher z-transformed Pearson correlation between predicted and actual BOLD');
ylabel('Fisher z-transformed Pearson correlation coefficient');


figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'theory', 'sprite', 'interaction', 'termination'});
cmap = [1 0.8 0.6 0.4]' * h(1).FaceColor;
plot_gp_CV_rois_helper(fs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names, [], cmap); %colormap(winter(3)));
title('Fisher z-transformed Pearson correlation between predicted and actual BOLD');
ylabel('Fisher z-transformed Pearson correlation coefficient');

figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'DQN', 'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
cmap = [1 0.9 0.8 0.7 0.6 0.5]' * h(2).FaceColor;
h = plot_gp_CV_rois_helper(fs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names, [], cmap); %colormap(autumn(5)));
title('Fisher z-transformed Pearson correlation between predicted and actual BOLD');
ylabel('Fisher z-transformed Pearson correlation coefficient');
%}


%{
figure;
plot_gp_CV_rois_helper(zs, 'ttest', regressor_names, roi_names);
title('Pearson correlation');
%}

%{
figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'sprite', 'interaction', 'termination'});
h = plot_gp_CV_rois_helper(zs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names);
%ylim([0 0.1]);
title('Fisher z-transformed Pearson correlation between predicted and actual BOLD');
ylabel('Fisher z-transformed Pearson correlation coefficient');
cmap = colormap(winter(3));
for i = 1:length(h)
    h(i).FaceColor = cmap(i,:);
end
%}



