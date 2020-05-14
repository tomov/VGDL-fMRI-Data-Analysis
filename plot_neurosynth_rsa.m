% plot permutation test results after neurosynth_rsa.m

%load mat/neurosynth_rsa_3_us=0_l=1_nperms=0_nroi=351.mat
%load mat/neurosynth_rsa_1_us=0_l=1_pi=355.mat
load mat/neurosynth_rsa_5_us=0_l=1_nperms=0_nroi=351.mat

ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T, all_subject_rhos, roi_masks);

names = region(which);

r = ceil(sqrt(length(roi_masks)));
c = ceil(length(roi_masks)/r);

figure;
for i = 1:length(roi_masks)
    subplot(r,c,i);

    hold on;
    hist(squeeze(perm_Rhos(i,1,:)));

    yl = ylim;
    line([Rho(i) Rho(i)], yl, 'color', 'red');
    text(Rho(i) * 1.1, yl(2) * 0.8, sprintf('p = %.3f', pval(i)));

    title(num2str(names(i)));
    xlabel('Spearman $\rho$', 'interpreter', 'latex');
end
