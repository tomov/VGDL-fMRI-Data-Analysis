% plot permutation test results after neurosynth_rsa.m

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
