% plot permutation test results after neurosynth_rsa_HRR.m


close all;
clear all;

load('mat/neurosynth_rsa_HRR.mat');

if contains(EXPT.rsadir, '_nosmooth')
    EXPT = vgdl_expt_nosmooth();
else
    EXPT = vgdl_expt();
end



alpha = 0.3; % signifinace level for thresholding (uncorr.) based on permutation tests
idx = find([ROI.p_rho] < alpha);


V = spm_vol('masks/mask.nii');

roi_masks = roi_masks(idx);
ROI = ROI(idx);

rho_map = nan(V.dim);
for r = 1:length(roi_masks)
    roi_mask = roi_masks{r};

    rho_map(roi_mask) = ROI(r).z_rho;
end
bspmview_wrapper(EXPT, rho_map);




figure;

subplot(1,3,1);

hold on;
hist(null_sym);
yl = ylim;
xl = xlim;
line([sym sym], yl, 'color', 'red');
text(xl(1), yl(2) * 0.9, sprintf('p = %.3f', p_sym));
set(gca, 'ytick', []);
title('Symmetry');
xlabel('symmetry coefficient', 'interpreter', 'latex');


subplot(1,3,2);
hist(null_diff);
yl = ylim;
line([diff diff], yl, 'color', 'red');
text(diff * 0.9, yl(2) * 0.9, sprintf('p = %.3f', p_diff));
set(gca, 'ytick', []);
title('Across - within game');
xlabel('$\Delta$ r', 'interpreter', 'latex');


subplot(1,3,3);
hist(null_rho);
yl = ylim;
line([rho rho], yl, 'color', 'red');
text(rho * 0.9, yl(2) * 0.9, sprintf('p = %.3f', p_rho));
set(gca, 'ytick', []);
title('HRR RSA match');
xlabel('Spearman $\rho$', 'interpreter', 'latex');
