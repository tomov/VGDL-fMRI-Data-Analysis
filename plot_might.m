% plot results from might.m

close all;
%clear all;

load('mat/might_knn_rsa=6_us=0_r=6.6667_znone.mat')
%load('mat/might_lda_shrinkage_rsa=5_us=0_r=2.6667_znone.mat');
%load('mat/might_knn_rsa=5_us=0_r=6.6667_zrun.mat');
%load('mat/might_knn_rsa=5_us=0_r=6.6667_znone.mat');
%load('mat/might_knn_rsa=5_us=1_r=6.6667_znone.mat');
%load('mat/might_knn_rsa=5_us=1_r=6.6667_znone.mat');
%load('mat/might_svm_linear_rsa=5_us=0_r=6.6667_znone.mat');
%load('mat/might_lda_shrinkage_rsa=1_us=0_r=6.6667_znone.mat');

% compute accuracy map
amap(mask) = mean(ams,1);

% count map (uncorr.)
cmap = zeros(size(mask));
%cmap(mask) = sum(pms < 0.05, 1); -- uncorr.

% count map (corr.)
alpha = 0.05;
howmany = zeros(nvoxels, 1); % for each voxel, how many subjects have it significant (according to pFDR)
for s = 1:length(subjects)
    [~, qm] = mafdr(pms(s,:)'); % Storey (2002)
    howmany = howmany + (qm < alpha);
end
cmap(mask) = howmany; % -- FDR

% t map
[h, p, ci, stats] = ttest(ams, 1/6);
tmap = zeros(size(mask));
tmap(mask) = stats.tstat;
tmap(tmap < 0) = 0; % o/w bspmview crashes

%
% bspmview boilerplate
%


if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

bspmview_wrapper(EXPT, tmap);

