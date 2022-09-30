clear all;
close all;


agg_filename = '/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_by_game_alpha=0.050_atlas=AAL2_GP_EMPA2_neuron.mat';
agg_filename

load(agg_filename);


game_names = {'Chase','Helper','Bait','Lemmings',{'Plaque', 'Attack'}, {'Avoid', 'George'},'Zelda'};
game_ids = {1, 2, 3, 4, 5, 5, 6};
subj_ids = {1:32, 1:32, 1:32, 1:32, 1:11, 12:32, 1:32};


%
%% fraction significant voxels
%
figure('pos', [99 201 1445 698]);
tiledlayout(length(game_names), 1, 'TileSpacing', 'none', 'Padding', 'none');

for g = 1:length(game_names)
    nexttile

    f = squeeze(fs(game_ids{g},:,:,subj_ids{g}));
    h = plot_gp_CV_rois_helper(f, 'signrank', 'median', regressor_names, roi_names, [], [], 10);
    title(game_names{g});
    ylabel('Fraction significant voxels');

end

