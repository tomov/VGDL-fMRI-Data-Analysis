% Insert the subject's SPM .mat filename here
EXPT = vgdl_expt();
model = 1;
subj = 1;

modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
spm_mat_file = fullfile(modeldir,'SPM.mat');

% Start batch
clear matlabbatch;
matlabbatch{1}.spm.util.voi.spmmat  = cellstr(spm_mat_file);
matlabbatch{1}.spm.util.voi.adjust  = 1;                    % Effects of interest contrast number
matlabbatch{1}.spm.util.voi.session = 1;                    % Session index
matlabbatch{1}.spm.util.voi.name    = 'test_voi';               % VOI name

matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii'};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.1;
%{
% Define thresholded SPM for finding the subject's local peak response
matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat      = {''};
matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast    = 1;     % Index of contrast for choosing voxels
matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc  = 'none';
matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh      = 0.001;
matlabbatch{1}.spm.util.voi.roi{1}.spm.extent      = 0;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask ...
    = struct('contrast', {}, 'thresh', {}, 'mtype', {});

% Define large fixed outer sphere
matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre     = [-15 -37 59]; % Set coordinates here
matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius     = 15;           % Radius (mm)
matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;

% Define smaller inner sphere which jumps to the peak of the outer sphere
matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre           = [0 0 0]; % Leave this at zero
matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius           = 8;       % Set radius here (mm)
matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm  = 1;       % Index of SPM within the batch
matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';    % Index of the outer sphere within the batch

% Include voxels in the thresholded SPM (i1) and the mobile inner sphere (i3)
matlabbatch{1}.spm.util.voi.expression = 'i1 & i3'; 
%}
matlabbatch{1}.spm.util.voi.expression = 'i1'; 

% optionally visualize result
matlabbatch{1}.spm.util.voi.disp = 0; 

% Run the batch
spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);
