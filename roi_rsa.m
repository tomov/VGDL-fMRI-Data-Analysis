function roi_rsa(rsa_idx, use_smooth, lateralized, nperms, roi_names, subbatch_size)

% RSA with custom ROIs
% copied from neurosynth_rsa.m

printcode;

if ~exist('subbatch_size', 'var')
    subbatch_size = 50; % don't do all ROIs at once; we OOM b/c all the Neural RDMs; too few though, and you end up loading betas too often
end

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

% get ROIs
[roi_masks, region] = get_apriori_rois(lateralized);
which = ismember(region, roi_names);
roi_masks = roi_masks(which);

tmp = join(cellfun(@num2str, region, 'UniformOutput', false), ',');
filename = sprintf('mat/roi_rsa_%d_us=%d_l=%d_nperms=%d_nroi=%d_rois=%s.mat', rsa_idx, use_smooth, lateralized, nperms, length(roi_masks), tmp{1});
filename

% run actual rsa
rsa_helper(EXPT, rsa_idx, roi_masks, filename, nperms, subbatch_size);
