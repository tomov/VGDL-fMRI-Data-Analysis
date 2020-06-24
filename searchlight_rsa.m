function searchlight_rsa(rsa_idx, use_smooth, nperms, subbatch_size, corr_type, dist, radius)

% RSA with ROIs from searchlights 
% copied from neurosynth_rsa.m

printcode;

%rsa_idx = 1;
%lateralized = true;
%use_smooth = true;
%nperms = 1000; % 0 for no permutation tests

if ~exist('subbatch_size', 'var')
    subbatch_size = 50; % don't do all ROIs at once; we OOM b/c all the Neural RDMs; too few though, and you end up loading betas too often
end

if ~exist('corr_type', 'var')
    corr_type = 'ktaub';
end

if ~exist('dist', 'var')
    dist = 'correlation';
end

if ~exist('radius', 'var')
    radius = 6 / 1.5; % mm -> voxels
else
    radius = radius / 1.5;
end

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

inds = 1:300000;

assert(nperms == 0);

filename = sprintf('mat/searchlight_rsa_%d_us=%d_nperms=%d_corr=%s_dist=%s_r=%.3f.mat', rsa_idx, use_smooth, nperms, corr_type, dist, radius);
filename

subjects = 1:length(EXPT.subject);

ccnl_rsa_searchlight(EXPT, rsa_idx, inds, radius, subbatch_size, subjects, corr_type, dist);

save(filename, '-v7.3');
