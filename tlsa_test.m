% This script demonstrates how to create and analyze a data set with TLSA
% momchil: c/p from might.m
% random effects 

addpath(genpath('../../tlsa_matlab'));

use_smooth = true;
rsa_idx = 3;

if use_smooth
    EXPT = vgdl_expt();
    maskfile = 'masks/mask.nii';
else
    EXPT = vgdl_expt_nosmooth();
    maskfile = 'masks/mask_nosmooth.nii';
end

clear all;
close all;

[mask, Vmask] = ccnl_load_mask(maskfile);

subjects = 1:length(EXPT.subject);

% Parameters of the synthetic data
N = 80;     % number of observations
C = 2;      % number of covariates
K = 10;     % number of latent sources
D = 3;      % number of features (3 spatial dimensions and 1 time dimension)
M = D+2;    % number of parameters per source
S = length(subjects);      % number of subjects
tau = 1;    % noise precision

% TLSA options (missing fields get set to defaults)
opts.mapfun = @(theta,R) map_st_rbf(theta,R);    % mapping function (spatiotemporal RBF)
opts.K = K;
opts.beta = 0.01;  % set to 0 to fit each subject independently

% Create the synthetic data set
%[r1 r2 r3 r4] = ndgrid(linspace(0,1,5)');
%R = [r1(:) r2(:) r3(:) r4(:)];  % location matrix
[x y z] = ind2sub(size(mask), find(mask));
R = cor2mni([x y z], Vmask.mat);

omega = randn(K,M);             % source parameters (here each subject has the same parameters)

rng(234);
vox = randsample(size(R,1), 5000); % subsample voxels, to make things faster

for s = 1:S
    subj = subjects(s);
    fprintf('    subj %d\n', subj);

    rsa = EXPT.create_rsa(rsa_idx, subj);

    modeldir = fullfile(EXPT.modeldir,['model',num2str(rsa.glmodel)],['subj',num2str(subj)]);
    spm_mat_file = fullfile(modeldir,'SPM.mat');
    load(spm_mat_file);

    assert(rsa_idx == 3);
    Y = ccnl_get_activations(EXPT, rsa.glmodel, mask, subj, true, true); % whiten & filter; see Diedrichsen et al. 2016
    Y = Y{1};
    %X = SPM.xX.xKXs.X; % whiten & high-pass filtered

    X = rsa.model(1).features;
    foldid = rsa.model(1).partitions;

    clear data;

    [~,truth] = max(X,[],2);
    pred = nan(size(truth));

    % CV
    folds = unique(foldid);
    for i = 1:length(folds)
        k = folds(i);

        % random effects for now (doesn't matter much; Gershman 2011)
        data(1).X = X(foldid ~= k, :);
        data(1).Y = Y(foldid ~= k, vox);
        data(1).R = R(vox,:);

        testdata(1).X = X(foldid == k, :);
        testdata(1).Y = Y(foldid == k, vox);
        testdata(1).R = R(vox,:);

        % run variational expectation-maximization algorithm

        results = tlsa_EM(data,opts);

        % decode covariates for test data
        %mu = tlsa_decode_gaussian(data,testdata,results);
        post = tlsa_decode_discrete(testdata,results, eye(6));
        post = post{1};

        [~, pred(foldid == k,:)] = max(post, [], 2);
    end

    acc(s) = mean(pred == truth);

end

% show inferred and ground truth covariates for a single subject
%figure;
%scatter(testdata(1).X(:),mu{1}(:)); lsline
%xlabel('Ground truth covariates','FontSize',15);
%ylabel('Decoded covariates','FontSize',15);
