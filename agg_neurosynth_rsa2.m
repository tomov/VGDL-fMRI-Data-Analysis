% aggregate results from neurosynth_rsa2.m
clear all;

use_smooth = true;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);
%subjects = [1 2 3 9 11 15 23 32];
model_name = 'EMPA';
what = 'theory';
glmodel = 9;
project = 0;
use_smooth = true;
maskname = 'mask';
neural_distance = 'correlation';

[mask_format, mask, Vmask] = get_mask_format_helper('masks/mask.nii'); % TODO 

agg_filename = sprintf('/Volumes/fMRI-2/Mac_mat/agg_neurosynth_rsa2_us=%d_glm=%d_model=%s_%s_project=%d.mat', use_smooth, glmodel, model_name, what, project);
agg_filename

for s = 1:length(subjects)
    subj_id = subjects(s);

    filename = sprintf('/Volumes/fMRI-2/Mac_mat/neurosynth_rsa2_subj=%d_us=%d_glm=%d_mask=%s_model=%s_%s_nsamples=100_project=%d_dist=%s.mat', subj_id, use_smooth, glmodel, maskname, model_name, what, project, neural_distance);
    filename
    load(filename, 'rho', 'roi_masks');

    if s == 1
        rhos = nan(length(subjects), sum(mask(:)));
    end

    for i = 1:length(roi_masks)
        roi_mask = roi_masks{i};
        assert(all(mask(roi_mask)));

        rhos(s, roi_mask(mask)) = rho(i);
    end
end

zs = atanh(rhos);

[h,p,ci,stats] = ttest(zs);
ts = stats.tstat;


tmap = zeros(size(mask));
tmap(mask) = ts;



agg_filename
save(agg_filename);


