function bspmview_wrapper(EXPT, map)

filename = bspmview_save_map(EXPT, map)

% view map
struc = fullfile(EXPT.modeldir,'mean.nii');
if exist(struc,'file')
    bspmview(filename, struc);
else
    disp('Using default base image');
    bspmview(filename);
end
