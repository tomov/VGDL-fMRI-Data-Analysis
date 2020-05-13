function [roi_masks, region] = get_apriori_rois(lateralized)

    % get neurosynth parcellation ROIs, normalized to subject mask
    %
    % returns:
    % roi_mask = cell array of 3D masks 
    % region = cell array of region names
    
    mask_filenames = { ...
        'masks/hippocampus.nii', ...
        'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii', ...
        };
    region = { ...
        'HC', ...
        'glm3', ...
        };

    group_mask = ccnl_load_mask('masks/mask.nii');
    
    for i = 1:length(mask_filenames)
        roi_masks{i} = ccnl_load_mask(mask_filenames{i});
        assert(all(size(roi_masks{i}) == size(group_mask)));
    end

