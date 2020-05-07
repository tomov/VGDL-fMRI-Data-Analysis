function [sphere_mask, sphere_coords] = create_spherical_mask_helper(mask, x, y, z, r, min_x, max_x, min_y, max_y, min_z, max_z, Vmask)

% does the bulk of the work in create_spherical_mask
% convenient to use if getting multiple masks in bulk (saves some overhead)
% Note that the coordinates are in native space, NOT in MNI space!
%
sphere_coords = [];

sphere_mask = zeros(size(mask));

for newx = floor(x - r) : ceil(x + r)
    if newx < min_x || newx > max_x, continue; end
    for newy = floor(y - r) : ceil(y + r)
        if newy < min_y || newy > max_y, continue; end
        for newz = floor(z - r) : ceil(z + r)
            if newz < min_z || newz > max_z, continue; end
            if ~mask(newx, newy, newz), continue; end
            if (x - newx)^2 + (y - newy)^2 + (z - newz)^2 > r^2, continue; end
            sphere_mask(newx, newy, newz) = 1;
           % mni = cor2mni([newx newy newz], Vmask.mat);
           % sphere_coords = [sphere_coords; mni];
            sphere_coords = [sphere_coords; newx newy newz];
        end
    end
end

sphere_mask = logical(sphere_mask);

