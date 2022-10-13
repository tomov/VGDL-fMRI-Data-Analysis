close all;
clear all;

EXPT = vgdl_expt();

[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');


glmodel = 102;
filename = sprintf('get_average_beta_glm=%d.mat', glmodel);
load(fullfile(get_mat_dir(2), filename), 'unique_regressors', 'avg_B');


ix = find(contains(unique_regressors, 'theory_change_flag'));

bmap = zeros(size(whole_brain_mask));
bmap(whole_brain_mask) = avg_B(ix,:);

bspmview_wrapper(EXPT, bmap);
