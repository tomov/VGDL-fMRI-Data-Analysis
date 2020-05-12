% functional connectivity

subj = 1;
seed_mask = 'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii';
EXPT = vgdl_expt();

seed = ccnl_get_activations(EXPT, 1, seed_mask, subj);
seed = seed{1};

seed = mean(seed, 2); % TODO eigenvariate

group_mask = 'masks/mask.nii';

disp('getting B');
tic
B = ccnl_get_activations(EXPT, 1, group_mask, subj);
B = B{1};
toc

seed = zscore(seed, 0, 1);
B = zscore(B, 0, 1);

disp('comupting correlations');
tic
% corr = dot product of z-scored vectors
r = seed' * B / size(B,1);
toc

z = atanh(r); % fisher z transform

[m, V, Y] = ccnl_load_mask('masks/mask.nii');
V.fname = 'masks/temp.nii'; % change immediately!
zmap = zeros(size(Y));
zmap(m) = z;

spm_write_vol(V, zmap);
struc = fullfile(EXPT.modeldir, 'mean.nii');

bspmview(V.fname, struc);
