function view(filename)

EXPT = vgdl_expt();
struc = fullfile(EXPT.modeldir,'mean.nii');
bspmview(filename, struc);

