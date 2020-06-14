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
%load('mat/neurosynth_rsa_6_us=1_l=1_nperms=1000_nroi=351.mat');

%load('mat/roi_rsa_6_us=1_l=1_nperms=1000_nroi=23_all.mat');

%load mat/neurosynth_rsa_5_us=0_l=1_nperms=0_nroi=351.mat

load mat/neurosynth_rsa_5_us=0_l=1_nperms=500_nroi=351.mat

%load mat/neurosynth_rsa_1_us=1_l=1_nperms=1000_nroi=351.mat

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
%tbl = table(Rho(i), T(i), names(i)');
tbl = table(Rho(i), T(i), names(i));

alpha = 0.1; % signifinace level for thresholding (uncorr.) based on permutation tests
idx = find(pval < alpha);

ccnl_rsa_view(EXPT, rsa_idx, 1, Rho_adj(idx), all_subject_rhos(idx,:,:), roi_masks(idx));


%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T, all_subject_rhos, roi_masks);
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, Rho_adj(pval < 1.1), all_subject_rhos(pval < 1.1,:,:), roi_masks(pval < 1.1));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T(pval < 0.1), all_subject_rhos(pval < 0.1,:,:), roi_masks(pval < 0.1));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, Rho_adj(pval < 0.1), all_subject_rhos(pval < 0.1,:,:), roi_masks(pval < 0.1));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T(T > -1.3), all_subject_rhos(T > -1.3), roi_masks(T > -1.3));
%ccnl_rsa_view(vgdl_expt(), rsa_idx, 1, T(T > -1.3), all_subject_rhos(T > -1.3), roi_masks(T > -1.3));


% compare p-values from different group-level (one-sided) hypothesis tests
%
assert(size(all_subject_rhos,2) == 1);
n = size(all_subject_rhos, 3);
F = atanh(all_subject_rhos); % fisher z score
z = sqrt((n - 3) / 1.06) * F; % z-statistic; https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient

p_perm = pval';
[~,p_t] = ttest(squeeze(F)', 0, 'tail', 'right');
[~,p_z] = ztest(squeeze(z)', 0, 1, 'tail', 'right');

figure;

subplot(1,3,1);
hold on;
scatter(p_perm, p_t);
plot([0 1], [0 1]);
xlabel('permutation test p-value');
ylabel('t test p-value');

subplot(1,3,2);
hold on;
scatter(p_perm, p_z);
plot([0 1], [0 1]);
xlabel('permutation test p-value');
ylabel('z test p-value');

title('each data point = neurosynth ROI; keep in mind those are dependent data points (b/c brain regions are dependent)');

subplot(1,3,3);
hold on;
scatter(p_z, p_t);
plot([0 1], [0 1]);
xlabel('z test p-value');
ylabel('t test p-value');



[r,p] = corr(p_perm', p_t');
fprintf('p_perm vs. p_t: r = %.3f, p = %e\n', r, p);

[r,p] = corr(p_perm', p_z');
fprintf('p_perm vs. p_z: r = %.3f, p = %e\n', r, p);

[r,p] = corr(p_z', p_t');
fprintf('p_z vs. p_t: r = %.3f, p = %e\n', r, p);


r = ceil(sqrt(length(idx)));
c = ceil(length(idx)/r);

%{
figure;
for i = 1:length(idx)
    subplot(r,c,i);

    hold on;
    hist(squeeze(perm_Rhos(idx(i),1,:)));

    yl = ylim;
    line([Rho(idx(i)) Rho(idx(i))], yl, 'color', 'red');
    text(Rho(idx(i)) * 1.1, yl(2) * 0.8, sprintf('p = %.3f', pval(idx(i))));
    set(gca, 'ytick', []);

    %title(names{idx(i)});
    title(num2str(names(idx(i))));
    xlabel('Spearman $\rho$', 'interpreter', 'latex');
end
%}
