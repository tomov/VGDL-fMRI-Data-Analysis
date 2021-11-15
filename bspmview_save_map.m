function filename = bspmview_save_map(EXPT, map, filename)

if ~exist('filename', 'var') || isempty(filename)
    filename = fullfile(EXPT.rsadir, ['temp_map.nii']);
end

%V = spm_vol(fullfile('masks', 'spmT_0001_pilot.nii'));
V = spm_vol(fullfile('masks', 'spmT_0001.nii')); 

% hacks to make it save the t-map as a t-map
V.fname = filename; % change immediately! SPM.mat file needs to be there too, to allow cluster FWE correction
V.dt = [16 0];
V.private.dat.dtype = 'FLOAT32-LE';
V.private.dat.fname = V.fname;

% save map
V.fname
spm_write_vol(V, map);


