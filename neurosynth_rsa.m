
% copied from Exploration / neurosynth_CV.m

printcode;

rsa_idx = 1;
lateralized = true;
use_smooth = false;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

% get ROIs
group_mask_filename = fullfile('masks', 'mask.nii');
if lateralized
    parcellation_file = fullfile('masks', 'Neurosynth Parcellation_2_lateralized.nii');
else
    parcellation_file = fullfile('masks', 'Neurosynth Parcellation_2.nii');
end

[~, Vparcel, parcellation_vol] = ccnl_load_mask(parcellation_file);
parcellation_vol = round(parcellation_vol);

if ~exist('parcel_idxs', 'var') || isempty(parcel_idxs)
    parcel_idxs = unique(parcellation_vol(:));
else
    assert(all(ismember(parcel_idxs, unique(parcellation_vol(:)))));
end


filename = sprintf('mat/neurosynth_rsa_%d_us=%d_l=%d_pi=%d.mat', rsa_idx, use_smooth, lateralized, length(parcel_idxs));
disp(filename);


region = [];
roi_masks = {};
for i = 1:length(parcel_idxs)

    parcel_idx = parcel_idxs(i);
    if parcel_idx == 0
        continue;
    end

    i
    region = [region; parcel_idx];
    fprintf('parcel = %d\n', parcel_idx);

    % ROI
    mask = parcellation_vol == parcel_idx;

    % normalize mask
    group_vol = spm_vol(group_mask_filename);
    group_mask = spm_read_vols(group_vol);

    [x, y, z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
    cor = mni2cor(cor2mni([x y z], Vparcel.mat), group_vol.mat); % voxel coords in AAL2 space --> voxel coords in MNI space --> voxel coords in our space
    ind = sub2ind(size(group_mask), cor(:,1), cor(:,2), cor(:,3)); % voxel coords in our space --> voxel indices
    
    % Reconstruct mask in our space
    %
    Vmask = group_vol;
    Vmask.fname = 'tmp.nii'; % CRUCIAL! o/w overwrite mask.nii
    mask = zeros(size(group_mask));
    mask(ind) = 1; % voxel indices --> binary mask
    
    % Only include voxels that are part of the subject group-level mask
    % i.e. that have non-NaN betas for all subjects
    %
    mask = mask & group_mask;

    if sum(mask(:)) < 5
        disp('skipping parcel -- small or empty mask');
        continue; % some masks are already lateralized
    end

    roi_masks = [roi_masks; {mask}];
end

save(filename, '-v7.3');

% run RSA
%
[Rho, H, T, P, all_subject_rhos, Behavioral, Neural] = ccnl_rsa(EXPT, rsa_idx, roi_masks);

save(filename, '-v7.3');

% view RSA results
%
ccnl_rsa_view(EXPT, rsa_idx, 1, T, all_subject_rhos, roi_masks);
