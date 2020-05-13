function view_rois(roi_masks)

    % helper f'n to plot rois real quick

    EXPT = vgdl_expt();
    [m, V, Y] = ccnl_load_mask('masks/mask.nii');
    V.fname = 'masks/temp.nii'; % change immediately!
    V.private.dat.dtype = 'FLOAT32-LE';
    V.dt = [16 0];

    map = zeros(size(Y));
    for i = 1:length(roi_masks)
        map(roi_masks{i}) = i;
    end

    spm_write_vol(V, map);
    struc = fullfile(EXPT.modeldir, 'mean.nii');

    bspmview(V.fname, struc);
