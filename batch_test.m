%-----------------------------------------------------------------------
% Job saved on 14-May-2020 12:24:26 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
function batch_test(EXPT, subj)

%EXPT = vgdl_expt();
model = 3;
%subj = 1;
run = 1;

modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
spm_mat_file = fullfile(modeldir,'SPM.mat');


% Start batch
% see matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/data/klinik/DOC07/crone/dcmmodels/masks/DLT.nii'};
% and https://en.wikibooks.org/wiki/SPM/Timeseries_extraction
%{
clear matlabbatch;
matlabbatch{1}.spm.util.voi.spmmat  = cellstr(spm_mat_file);
matlabbatch{1}.spm.util.voi.adjust  = 0;                    % Effects of interest contrast number
matlabbatch{1}.spm.util.voi.session = run;                    % Session index
matlabbatch{1}.spm.util.voi.name    = 'test_voi';               % VOI name

matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii'};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.1;

matlabbatch{1}.spm.util.voi.expression = 'i1'; 

% optionally visualize result
matlabbatch{1}.spm.util.voi.disp = 0; 

% Run the batch
spm('defaults', 'FMRI');
voi_job = spm_jobman('run',matlabbatch);
%}

job.spmmat  = cellstr(spm_mat_file);
job.adjust  = 0;                    % Effects of interest contrast number
job.session = run;                    % Session index
job.name    = 'test_voi';               % VOI name

job.roi{1}.mask.image = {'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii'};
job.roi{1}.mask.threshold = 0.1;

job.expression = 'i1'; 

out = spm_run_voi(job);
%load(voi_job{1}.voimat{1}, 'xY');
load(out.voimat{1});
load(spm_mat_file);

ppi = spm_peb_ppi(SPM, 'sd', xY, [], 'ppibro', true);


%{
clear matlabbatch;

matlabbatch{1}.spm.stats.ppi.spmmat = {spm_mat_file};

%matlabbatch{1}.spm.stats.ppi.type.sd.voi = {fullfile(modeldir, 'VOI_test_voi_1.mat')};
matlabbatch{1}.spm.stats.ppi.type.sd.voi = {xY};

%matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {'/Users/momchil/Dropbox/Research/demo/attention/GLM/VOI_bvoi_name_1.mat'};
%matlabbatch{1}.spm.stats.ppi.type.ppi.u = [2 1 -1
%                                           3 1 1];

matlabbatch{1}.spm.stats.ppi.name = 'test_sd';
matlabbatch{1}.spm.stats.ppi.disp = 1;

spm_jobman('run',matlabbatch);
%}
