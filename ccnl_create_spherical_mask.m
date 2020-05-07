function [sphere_mask, sphere_coords, sphere_vol] = ccnl_create_spherical_mask(x, y, z, r, filename, intersect_mask)

% Create a spherical mask around a voxel. 
%
% USAGE:
%   [sphere_mask, sphere_coords] = ccnl_create_spherical_mask(x, y, z, r, filename, intersect_mask)
%
% INPUT:
%   x, y, z = coordinates of the center voxel in native coordinate space 
%             (use mni2cor first if they're in MNI space)
%   r = radius in voxels i.e. the native coordinate space, NOT MNI space!
%   filename = (optional) output filename where to save the .nii file. If not provided,
%              file won't be saved
%   intersect_mask = (optional) intersect the sphere with given mask (e.g. cluster from contrast) 
%
% OUTPUT:
%   sphere_mask = the resulting mask as a 3D binary vector
%   sphere_coords = V x 3 list of MNI coordinates of the sphere voxels   
%   sphere_vol = the corresponding SPM volume variable
%

% load whole-brain mask so we don't pick voxels outside of it
%
[mask, Vmask] = ccnl_load_mask('masks/mask.nii');
Vmask.fname = ''; % !!!! IMPORTANT!!!!! b/c we return it
if exist('filename', 'var') && ~isempty(filename)
    Vmask.fname = filename;
end
sphere_vol = Vmask;

% convert sphere center to native coordinate space
%
%cor = mni2cor([x y z], Vmask.mat);
%x = cor(1);
%y = cor(2);
%z = cor(3);

% find boundary coordinates
%
[all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
min_x = min(all_x);
max_x = max(all_x);
min_y = min(all_y);
max_y = max(all_y);
min_z = min(all_z);
max_z = max(all_z);


% create the spherical mask
%
[sphere_mask, sphere_coords] = create_spherical_mask_helper(mask, x, y, z, r, min_x, max_x, min_y, max_y, min_z, max_z, Vmask);

% optionally intersect with given mask
%
if exist('intersect_mask', 'var') && ~isempty(intersect_mask)
    sphere_mask = sphere_mask & intersect_mask;
end

% optionally save the mask
%
if exist('filename', 'var') && ~isempty(filename)
    sphere_vol.fname = filename;
    spm_write_vol(sphere_vol, sphere_mask);
end

