% Script to create masks for subcortical RL ROIs
%

% The ROI RL atlas labels (bilateral by default)
%
rois = {'Ca', 'Pu', 'NAC', 'GPe', 'GPi', 'SNr', 'STH', 'SNc', 'VTA'};

group_mask_filename = fullfile('masks', 'mask.nii');

% Create the masks
%
for roi = rois
    roi = roi{1};
    
    % Create and save the mask
    %
    ccnl_create_mask({roi}, fullfile('masks', [roi, '.nii']), 'RL', true, group_mask_filename);
    ccnl_create_mask({roi}, fullfile('masks', [roi, '_unnormalized.nii']), 'RL', false);
end
