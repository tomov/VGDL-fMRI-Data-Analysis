function masks = create_masks_from_mni(file_prefix, mni, do_create, sphere)

    % create spherical masks around MNI coordinates
    % 
    % file_prefix: prefix of the mask filenames where to save them
    % mni: [n x 3] list of MNI coordinates
    % do_create: whether to actually create the files (or just generate the filenames)
    % sphere: sphere radius in mm

    r = sphere / 1.5; % mm -> voxels

    [mask, V, Y] = ccnl_load_mask('masks/mask.nii');
    V.fname = 'masks/temp.nii'; % change immediately!

    for i = 1:size(mni, 1)
        cor = mni2cor(mni(i,:), V.mat);

        V.fname = fullfile('masks', sprintf('%s_%d_%d_%d_r=%.1fmm.nii', file_prefix, mni(i,1), mni(i,2), mni(i,3), sphere));

        if do_create
            [sphere_mask] = ccnl_create_spherical_mask(cor(1), cor(2), cor(3), r);
            sphere_mask = sphere_mask & mask; % exclude voxels outside of brain
            spm_write_vol(V, sphere_mask);
        end

        masks{i} = V.fname;
    end

