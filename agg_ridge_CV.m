% aggregate results from fit_ridge_CV.m
% copy of agg_gp_CV.m

clear all;

use_smooth = true;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);
model_name = 'EMPA';
what = 'theory';
glmodel = 1;
project = 1;
subsample = 0;
fasse_ncf = false;

[mask_format, mask, Vmask] = get_mask_format_helper('masks/mask.nii'); % TODO 

%agg_filename = sprintf('mat/agg_ridge_CV_us=%d_glm=%d_%s.mat', use_smooth, glmodel, what);
agg_filename = fullfile(get_mat_dir(), sprintf('agg_ridge_CV_us=%d_glm=%d_model=%s_%s_subsample=%d_project=%d.mat', use_smooth, glmodel, model_name, what, subsample, project));
agg_filename

for s = 1:length(subjects)
    subj_id = subjects(s);

    %filename = sprintf('mat/fit_ridge_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_%s.mat', subj_id, use_smooth, glmodel, what);
    %filename = sprintf('mat/fit_ridge_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_subsample=1_%s.mat', subj_id, use_smooth, glmodel, what);
    filename = fullfile(get_mat_dir(fasse_ncf), sprintf('fit_ridge_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_subsample=%d_project=%d_saveYhat=0.mat', subj_id, use_smooth, glmodel, model_name, what, subsample, project));
    filename
    load(filename);

    % MSEs
    MSEs(s,:) = mean(MSE_CV, 1);
    SMSEs(s,:) = mean(SMSE_CV, 1);

    % R2s
    R2s(s,:) = mean(R2_CV, 1);
    adjR2s(s,:) = mean(adjR2_CV, 1);;

    % CV Pearson r's across subjects
    rs(s,:) = mean(r_CV, 1);
end

% Fisher z transform
zs = atanh(rs);

% t-test Pearson corr across subjects
[h,p,ci,stats] = ttest(zs);
%[h,p,ci,stats] = ttest(zs, null_zs);
ts = stats.tstat;

% create SPMs
%
R2map = zeros(size(mask));
R2map(mask) = mean(R2s,1);

adjR2map = zeros(size(mask));
adjR2map(mask) = mean(adjR2s,1);

tmap = zeros(size(mask));
tmap(mask) = ts;



agg_filename
save(agg_filename);

