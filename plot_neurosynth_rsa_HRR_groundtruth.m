% plot permutation test results after neurosynth_rsa_HRR_groundtruth.m


close all;
clear all;

%load('mat/neurosynth_rsa_HRR_groundtruth_us=1_nperm=1000.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_us=0_nperm=1000.mat');
%load('mat/neurosynth_rsa_HRR_groundtruth_us=0_nperm=0.mat');
load('mat/neurosynth_rsa_HRR_groundtruth_us=0_nperm=1000_game.mat');

if contains(EXPT.rsadir, '_nosmooth')
    EXPT = vgdl_expt_nosmooth();
else
    EXPT = vgdl_expt();
end

what = 'rho'; % options: 'rho', 'sym', 'diff'


%% spearman vs. null distribution group-level stats
% sanity check, to make sure permuting wasn't screwed up
% works for what = 'rho' only
%{

p_t = [ROI.(['p_', what, '_t'])];
p_z = [ROI.(['p_', what])];

figure;
hold on;
scatter(p_z, p_t);
plot([0 1], [0 1]);
hold off;
xlabel('z test p-value');
ylabel('t test p-value');
%}



%% brain map


alpha = 0.10; % signifinace level for thresholding (uncorr.) based on permutation tests
idx = find([ROI.(['p_', what])] < alpha);
%idx = find([ROI.(['p_', what, '_t'])] < alpha);

V = spm_vol('masks/mask.nii');

roi_masks = roi_masks(idx);
ROI = ROI(idx);

map = nan(V.dim);
for r = 1:length(roi_masks)
    roi_mask = roi_masks{r};

    map(roi_mask) = ROI(r).(['z_', what]);
    %map(roi_mask) = ROI(r).([what, '_tstat']).tstat;
end
bspmview_wrapper(EXPT, map);





%% null distributions histograms


null_sym = ROI(1).null_sym;
null_diff = ROI(1).null_diff;
null_rho = ROI(1).null_rho;
sym = ROI(1).sym;
diff = ROI(1).diff;
rho = ROI(1).rho;
p_sym = ROI(1).p_sym;
p_diff = ROI(1).p_diff;
p_rho = ROI(1).p_rho;


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
