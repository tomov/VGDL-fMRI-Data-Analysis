function neurosynth_rsa(rsa_idx, use_smooth, lateralized, nperms, parcel_idx, subbatch_size)

% RSA with ROIs from neurosynth parcellation
% copied from Exploration / neurosynth_CV.m

printcode;

%rsa_idx = 1;
%lateralized = true;
%use_smooth = true;
%nperms = 1000; % 0 for no permutation tests

if ~exist('subbatch_size', 'var')
    subbatch_size = 50; % don't do all ROIs at once; we OOM b/c all the Neural RDMs; too few though, and you end up loading betas too often
end

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

% get ROIs
[roi_masks, region] = get_neurosynth_rois(lateralized);

if ~exist('parcel_idx', 'var') || isempty(parcel_idx)
    filename = sprintf('mat/neurosynth_rsa_%d_us=%d_l=%d_nperms=%d_nroi=%d.mat', rsa_idx, use_smooth, lateralized, nperms, length(roi_masks));
else
    which = ismember(region, parcel_idx);
    roi_masks = roi_masks(which);

    tmp = join(cellfun(@num2str, num2cell(parcel_idx), 'UniformOutput', false), ',');
    filename = sprintf('mat/neurosynth_rsa_%d_us=%d_l=%d_nperms=%d_nroi=%d_pi=%s.mat', rsa_idx, use_smooth, lateralized, nperms, length(roi_masks), tmp{1});
end
filename

% run actual rsa
rsa_helper(EXPT, rsa_idx, roi_masks, filename, nperms, subbatch_size);
