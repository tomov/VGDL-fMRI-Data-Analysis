% plot permutation test results after neurosynth_rsa.m

close all;
clear all;

%load mat/neurosynth_rsa_3_us=0_l=1_nperms=0_nroi=351.mat
%load mat/neurosynth_rsa_1_us=0_l=1_pi=355.mat
%load mat/neurosynth_rsa_5_us=0_l=1_nperms=0_nroi=351.mat
%load mat/neurosynth_rsa_6_us=1_l=1_nperms=0_nroi=351.mat

%load mat/neurosynth_rsa_1_us=0_l=1_nperms=1000_nroi=351.mat
%load mat/neurosynth_rsa_1_us=1_l=1_nperms=1000_nroi=351.mat
%load mat/neurosynth_rsa_6_us=0_l=1_nperms=0_nroi=351.mat

%load('mat/neurosynth_rsa_5_us=0_l=1_nperms=200_nroi=18_pi=277,351,239,383,234,106,77,126,335,72,347,55,300,313,40,168,298,175.mat');
%load('mat/neurosynth_rsa_5_us=0_l=1_nperms=200_nroi=14_pi=149,95,207,56,45,91,22,264,274,112,94,305,169,209.mat');

%load('mat/neurosynth_rsa_6_us=0_l=1_nperms=1000_nroi=351.mat');
load('mat/neurosynth_rsa_6_us=1_l=1_nperms=1000_nroi=351.mat');

%load mat/neurosynth_rsa_5_us=0_l=1_nperms=0_nroi=351.mat

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

if contains(EXPT.rsadir, '_nosmooth')
    EXPT = vgdl_expt_nosmooth();
else
    EXPT = vgdl_expt();
end

if exist('region', 'var')
    % neurosynth
    if ~exist('which', 'var')
        which = logical(ones(size(region)));
    end
    names = region(which);
else
    % roi

end

[~,i] = sort(T);
tbl = table(Rho(i), T(i), names(i));

alpha = 0.01;
idx = find(pval < alpha);

%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T, all_subject_rhos, roi_masks);
ccnl_rsa_view(EXPT, rsa_idx, 1, Rho_adj(idx), all_subject_rhos(idx,:,:), roi_masks(idx));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, Rho_adj(pval < 1.1), all_subject_rhos(pval < 1.1,:,:), roi_masks(pval < 1.1));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T(pval < 0.1), all_subject_rhos(pval < 0.1,:,:), roi_masks(pval < 0.1));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, Rho_adj(pval < 0.1), all_subject_rhos(pval < 0.1,:,:), roi_masks(pval < 0.1));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T(T > -1.3), all_subject_rhos(T > -1.3), roi_masks(T > -1.3));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T(T > -1.3), all_subject_rhos(T > -1.3), roi_masks(T > -1.3));

r = ceil(sqrt(length(idx)));
c = ceil(length(idx)/r);

figure;
for i = 1:length(idx)
    subplot(r,c,i);

    hold on;
    hist(squeeze(perm_Rhos(idx(i),1,:)));

    yl = ylim;
    line([Rho(idx(i)) Rho(idx(i))], yl, 'color', 'red');
    text(Rho(idx(i)) * 1.1, yl(2) * 0.8, sprintf('p = %.3f', pval(idx(i))));
    set(gca, 'ytick', []);

    title(num2str(names(idx(i))));
    %title(num2str(i));
    xlabel('Spearman $\rho$', 'interpreter', 'latex');
end
