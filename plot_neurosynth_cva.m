% plot permutation test results after neurosynth_cva.m

close all;
clear all;

%load('mat/neurosynth_cva_1_us=0_l=1_cons=vgfmri3_bait+vgfmri3_chase+vgfmri3_helper+vgfmri3_lemmings+vgfmri3_plaqueAttack+vgfmri3_zelda_nroi=351.mat');
load('mat/neurosynth_cva_1_us=1_l=1_cons=vgfmri3_bait+vgfmri3_chase+vgfmri3_helper+vgfmri3_lemmings+vgfmri3_plaqueAttack+vgfmri3_zelda_nroi=351.mat');

if exist('region', 'var')
    % neurosynth
    if ~exist('which', 'var')
        which = logical(ones(size(region)));
    end
    names = region(which);
else
    % roi

end

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

mmin = 0;
idx = find(m > mmin);

% c/p ccnl_rsa_view

V = spm_vol('masks/mask.nii');

% hacks to make it save the t-map as a t-map
V.fname = fullfile(EXPT.rsadir, ['temp_cmap.nii']); % change immediately!
V.dt = [16 0];
V.private.dat.dtype = 'FLOAT32-LE';
V.private.dat.fname = V.fname;

cmap = nan(V.dim);
for i = 1:length(idx)
    roi_mask = roi_masks{idx(i)};

    % load roi mask
    [roi_mask_format, roi_mask] = get_mask_format_helper(roi_mask);
    assert(strcmp(roi_mask_format, 'mask'), 'Improper mask');

    cmap(roi_mask) = m(idx(i));
end

% save cmap
V.fname
spm_write_vol(V, cmap);

% view cmap
struc = fullfile(EXPT.modeldir,'mean.nii');
if exist(struc,'file')
    bspmview(V.fname, struc);
else
    bspmview(V.fname);
end
