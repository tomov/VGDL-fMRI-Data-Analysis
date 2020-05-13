function [roi_masks, region] = get_apriori_rois(lateralized)

    assert(lateralized);

    % get neurosynth parcellation ROIs, normalized to subject mask
    %
    % returns:
    % roi_mask = cell array of 3D masks 
    % region = cell array of region names
    
    mask_filenames = { ...
        'masks/v1_L.nii', ... 
        'masks/v1_R.nii', ... 
        'masks/m1_L.nii', ... 
        'masks/m1_R.nii', ... 
        'masks/s1_L.nii', ... 
        'masks/s1_R.nii', ... 
        'masks/hippocampus_L.nii', ... 
        'masks/hippocampus_R.nii', ... 
        'masks/mOFC_L.nii', ... 
        'masks/mOFC_R.nii', ... 
        'masks/lOFC_L.nii', ... 
        'masks/lOFC_R.nii', ... 
        'masks/aOFC_L.nii', ... 
        'masks/aOFC_R.nii', ... 
        'masks/pOFC_L.nii', ... 
        'masks/pOFC_R.nii', ... 
        'masks/rectus_L.nii', ... 
        'masks/rectus_R.nii', ... 
        'masks/vmPFC_L.nii', ... 
        'masks/vmPFC_R.nii', ... 
        'masks/IFGorb_L.nii', ... 
        'masks/IFGorb_R.nii', ... 
        'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii', ...
        };

    region = { ...
        'hippocampus_L', ... 
        'hippocampus_R', ... 
        'mOFC_L', ... 
        'mOFC_R', ... 
        'lOFC_L', ... 
        'lOFC_R', ... 
        'aOFC_L', ... 
        'aOFC_R', ... 
        'pOFC_L', ... 
        'pOFC_R', ... 
        'rectus_L', ... 
        'rectus_R', ... 
        'vmPFC_L', ... 
        'vmPFC_R', ... 
        'IFGorb_L', ... 
        'IFGorb_R', ...
        'glm3', ...
        };

    group_mask = ccnl_load_mask('masks/mask.nii');
    
    for i = 1:length(mask_filenames)
        roi_masks{i} = ccnl_load_mask(mask_filenames{i});
        assert(all(size(roi_masks{i}) == size(group_mask)));
    end

