% plot permutation test results after neurosynth_rsa.m

%load mat/neurosynth_rsa_3_us=0_l=1_nperms=0_nroi=351.mat
%load mat/neurosynth_rsa_1_us=0_l=1_pi=355.mat
load mat/neurosynth_rsa_5_us=0_l=1_nperms=0_nroi=351.mat

%load mat/roi_rsa_5_us=0_l=1_nperms=1000_nroi=23_all.mat
%load mat/roi_rsa_5_us=0_l=1_nperms=1000_nroi=23_all.mat
%load mat/roi_rsa_3_us=0_l=1_nperms=100_nroi=23_all.mat

%load('mat/neurosynth_rsa_1_us=1_l=1_nperms=10000_nroi=8_pi=177,294,365,38,293,174,194,118.mat')

%[roi_masks, region] = get_neurosynth_rois(true);
%parcel_idx = [177,146,293,190,47,341,140,230,353];
%which = ismember(region, parcel_idx);
%load('mat/neurosynth_rsa_3_us=0_l=1_nperms=200_nroi=9_pi=177,146,293,190,47,341,140,230,353.mat');

%load mat/neurosynth_rsa_3_us=0_l=1_nperms=0_nroi=351.mat
%load('mat/neurosynth_rsa_3_us=0_l=1_nperms=200_nroi=5_pi=99,77,137,51,341.mat');

[~,i] = sort(T);
tbl = table(Rho(i), T(i), region(i));

ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T, all_subject_rhos, roi_masks);
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T(T > -1.3), all_subject_rhos(T > -1.3), roi_masks(T > -1.3));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T(T > -1.3), all_subject_rhos(T > -1.3), roi_masks(T > -1.3));

if exist('region', 'var')
    % neurosynth
    names = region(which);
else
    % roi

end

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
    %title(num2str(i));
    xlabel('Spearman $\rho$', 'interpreter', 'latex');
end
