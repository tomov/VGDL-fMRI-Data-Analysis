function [voxelsToNeighbors, numberOfNeighbors] = might_computeNeighborsSphere(coords, r, mask)
% Same as the searchmight computeNeighboursWithinRadius except for spherical searchlight.
% Avoids voxels not in mask. Also avoids voxels that are NaN
% USAGE:
% [meta.voxelsToNeighbours, meta.numberOfNeighbours] = might_computeNeighborsSphere(meta.colToCoord, r, mask);

voxelsToNeighbors = zeros(size(coords,1), ceil(2*r+1)^3);
numberOfNeighbors = zeros(size(coords,1), 1);

% get mapping x,y,z --> idx
assert(sum(mask(:)) == size(coords, 1));
idxs = nan(size(mask)); % for reverse indexing
idxs(mask) = 1:sum(mask(:));

% sanity check mask voxel indexing
[x y z] = ind2sub(size(mask), find(mask));
assert(immse([x y z], coords) == 0);

% get boundaries
[all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
min_x = min(all_x);
max_x = max(all_x);
min_y = min(all_y);
max_y = max(all_y);
min_z = min(all_z);
max_z = max(all_z);

max_n = 0;

% for each voxel = sphere center
%
for i = 1:size(coords, 1)

    x = coords(i,1);
    y = coords(i,2);
    z = coords(i,3);

    neighbors = [];

    for newx = floor(x - r) : ceil(x + r)
        if newx < min_x || newx > max_x, continue; end
        for newy = floor(y - r) : ceil(y + r)
            if newy < min_y || newy > max_y, continue; end
            for newz = floor(z - r) : ceil(z + r)
                if newz < min_z || newz > max_z, continue; end
                if ~mask(newx, newy, newz), continue; end
                if (x - newx)^2 + (y - newy)^2 + (z - newz)^2 > r^2, continue; end

                idx = idxs(newx,newy,newz);
                assert(~isnan(idx));

                neighbors = [neighbors, idx];
            end
        end
    end

    n = numel(neighbors);
    voxelsToNeighbors(i,1:n) = neighbors;
    numberOfNeighbors(i) = n;
    max_n = max(max_n, n);
end

voxelsToNeighbors = voxelsToNeighbors(:,1:max_n);
