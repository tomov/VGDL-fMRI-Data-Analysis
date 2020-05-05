
% copied from Exploration / neurosynth_CV.m

clear all;

printcode;

rsa_idx = 1;
subbatch_size = 50; % don't do all ROIs at once; we OOM b/c all the Neural RDMs
lateralized = true;
use_smooth = true;
nperms = 1000; % 0 for no permutation tests

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

atlas.id = []; % create bspmview atlas
atlas.label = {};

region = [];
roi_masks = {};
for i = 1:length(parcel_idxs)

    parcel_idx = parcel_idxs(i);
    if parcel_idx == 0
        continue;
    end

    atlas.id = [atlas.id parcel_idx];
    atlas.label = [atlas.label {num2str(parcel_idx)}];

    i
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
    region = [region; parcel_idx];
end

%roi_masks = roi_masks(161); % zoom in on 1 ROI for permutation test

save(filename, '-v7.3');

% run RSA in (sub)batches
%
for s = 1:subbatch_size:length(roi_masks)
    e = min(length(roi_masks), s + subbatch_size - 1);
    fprintf('s:e = %d:%d\n', s, e);

    [Rho(s:e,:), H(s:e,:), T(s:e,:), P(s:e,:), all_subject_rhos(s:e,:,:)] = ccnl_rsa(EXPT, rsa_idx, roi_masks(s:e));
end

Behavioral = ccnl_behavioral_rdms(EXPT, rsa_idx); % for plotting

save(filename, '-v7.3');

% permutations
%
if nperms > 0
    rng('shuffle');

    create_rsa = EXPT.create_rsa;

    perm_Rhos = nan(size(Rho,1), size(Rho,2), nperms);
    for i = 1:nperms
        seed = randi(1000000); % notice this will get affected by the rng calls inside create_rsa, but that's okay
        fprintf('perm = %d\n', i);
        
        EXPT.create_rsa = @(rsa_idx, subj_id) create_rsa(rsa_idx, subj_id, seed); % overwrite create_rsa with one that randomly shuffles 
        % run permuted RSA in (sub)batches
        %
        for s = 1:subbatch_size:length(roi_masks)
            e = min(length(roi_masks), s + subbatch_size - 1);
            fprintf('        perm s:e = %d:%d\n', s, e);

            [perm_Rhos(s:e,:,i)] = ccnl_rsa(EXPT, rsa_idx, roi_masks(s:e));
        end
    end

    % p-value = P(Rho >= rho | null)
    pval = mean(perm_Rhos >= Rho, 3);

    % adjusted rho = subtract mean & divide by stdev (for plotting)
    Rho_adj = (Rho - mean(perm_Rhos, 3)) ./ std(perm_Rhos, 3);

    save(filename, '-v7.3');
end

% view RSA results
%
ccnl_rsa_view(EXPT, rsa_idx, 1, T, all_subject_rhos, roi_masks);
