% functional connectivity

glm = 3;
seed_mask = 'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii';
EXPT = vgdl_expt();
group_mask = 'masks/mask.nii';


[m, V, Y] = ccnl_load_mask('masks/mask.nii');
V.fname = 'masks/temp.nii'; % change immediately!
V.private.dat.dtype = 'FLOAT32-LE';
V.dt = [16 0];

nsubj = 8;

for subj = 1:nsubj

    seed = ccnl_get_residuals(EXPT, glm, seed_mask, subj);
    %seed = ccnl_get_activations(EXPT, 1, seed_mask, subj); % too correlated w/ white matter & all; way too many confounds
    seed = seed{1};

    seed = mean(seed, 2); % TODO eigenvariate

    disp('getting B');
    tic
    B = ccnl_get_residuals(EXPT, glm, group_mask, subj);
    %B = ccnl_get_activations(EXPT, 1, group_mask, subj);
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

    if ~exist('zs', 'var')
        zs = zeros(nsubj, length(z));
        zs(subj,:) = z;
    end
    %zmap = zeros(size(Y));
    %zmap(m) = z;

end

[h,p,ci,stat] = ttest(zs);

tmap = zeros(size(Y));
tmap(m) = stat.tstat;

spm_write_vol(V, tmap);
struc = fullfile(EXPT.modeldir, 'mean.nii');

bspmview(V.fname, struc);
