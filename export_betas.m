% export betas for MVPA analysis in python

rsa_idxs = [1 6 5];

use_smooth = false;

if use_smooth
    EXPT = vgdl_expt();
    maskfile = 'masks/mask.nii';
    suffix = 'smooth';
else
    EXPT = vgdl_expt_nosmooth();
    maskfile = 'masks/mask_nosmooth.nii';
    suffix = 'nosmooth';
end

subjects = 1:length(EXPT.subject);

[mask, Vmask] = ccnl_load_mask(maskfile);

for rsa_idx = rsa_idxs
    rsa_idx

    for s = 1:length(subjects)
        subj = subjects(s);
        fprintf('    subj %d\n', subj);

        rsa = EXPT.create_rsa(rsa_idx, subj);

        [B, names] = ccnl_get_beta_series(EXPT, rsa.glmodel, subj, rsa.event, mask);

        filename = sprintf('mat/beta_series_glm%d_subj%d_%s.mat', rsa.glmodel, subj, suffix);
        filename
        save(filename, 'B', 'names', 'mask', 'Vmask', '-v7.3');
    end
end
