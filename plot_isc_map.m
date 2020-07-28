% plot ISC maps from Daphne
% https://github.com/tomov/VGDL-fMRI-Python-Data-Analysis/blob/master/Multivariate_analyses/Compare_isc_analyses.ipynb

[mask, Vmask] = ccnl_load_mask('masks/mask.nii');
EXPT = vgdl_expt();

%load('mat/collapsed_corrs_blocks.mat');
%load('isc_analysis_mat/iscs_r_games.mat');
%load('isc_analysis_mat/iscs_r_blocks.mat');
load('isc_analysis_mat/iscs_r_levels.mat');

[h,p,ci,stat] = ttest(atanh(mydata));

map = zeros(size(mask));
map(mask) = stat.tstat;

bspmview_wrapper(EXPT, map);


