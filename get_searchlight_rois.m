function [ROI] = get_searchlight_rois(mask, Vmask, radius)

    % radius in mm, Vmask is spm volume info
    %

    [all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in native space

    ROI(length(all_x)).cor = 'dummy'; % to initialize the array

    for i = 1:length(all_x) % for every voxel
        x = all_x(i);
        y = all_y(i);
        z = all_z(i);

        if mod(i, 10000) == 0
            disp(i);
        end

        ROI(i).cor = [x y z];
        ROI(i).mni = cor2mni([x y z], Vmask.mat);
        ROI(i).voxel_idx = [];
        %ROI(i).voxel_cor = []; % too big

        for newx = floor(x - radius) : ceil(x + radius)
            if newx < 1 || newx > size(mask, 1)
                continue;
            end

            for newy = floor(y - radius) : ceil(y + radius)
                if newy < 1 || newy > size(mask, 2)
                    continue; 
                end

                for newz = floor(z - radius) : ceil(z + radius)
                    if newz < 1 || newz > size(mask, 3)
                        continue; 
                    end

                    if ~mask(newx, newy, newz)
                        continue; 
                    end

                    if (x - newx)^2 + (y - newy)^2 + (z - newz)^2 > radius^2
                        continue;
                    end

                    %ROI(i).voxel_cor = [ROI(i).voxel_cor; newx newy newz];
                    ROI(i).voxel_idx = [ROI(i).voxel_idx, sub2ind(size(mask), newx, newy, newz)];
                end
            end
        end

    end

