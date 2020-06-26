
%{
clear all;
close all;
%}

EXPT = vgdl_expt();
glmodel = 21; % control GLM
mask = 'masks/ROI_x=48_y=12_z=30_85voxels_Sphere6.nii';
subjects = 1:length(EXPT.subject);


% load mask
[mask_format, mask, Vmask] = get_mask_format_helper(mask);

addpath(genpath('/ncf/gershman/Lab/scripts/gpml'));


subj = subjects(1); % TODO for each subject

[Y, K, W, R] = load_subject(EXPT, glmodel, subj, mask, Vmask);

%X = rand(size(Y,1), round(size(Y,1) * 1.2));
X = [];
for i = 0.01:0.01:1
    X = [X; cos([1:1000]*i)'];
end
ker = X*X'; % temporary
ker = 

Y = R*K*W*Y;
ker = R*K*W*ker*W'*K'*R';

%ker = nearestSPD(ker);

y = Y(:,1); % TODO for each voxel

options = optimoptions('fmincon','SpecifyObjectiveGradient',true)
% TODO optimize -- compute invKi only once
%sigma_hat = fmincon(@(sigma) gp_loglik(ker, y, sigma, true), 1, -1, 0, [], [], [], [], [], options);


assert(size(ker,1) == size(Y,1));
n = size(ker,1);
L = cholcov(ker); % not chol b/c sometimes positive semi-definite
covhyp = L(triu(true(n))); % see covDiscrete.m

% GP hyperparams
hyp = struct('mean', zeros(1, n), 'cov', covhyp, 'lik', exp(sigma));

% GP
x = [1:size(Y,1)]';
[nlz, ~] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);


function [Y, K, W, R] = load_subject(EXPT, glmodel, subj, mask, Vmask)
    % load subject data
    % Y = raw BOLD
    % K = filter matrix
    % W = whitening matrix
    % R = residual forming matrix
    %

    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
    load(fullfile(modeldir,'SPM.mat'));
    assert(isempty(Vmask) || isequal(SPM.xY.VY(1).dim, Vmask.dim), 'Different dimensions between mask and activations');
    assert(ndims(mask) < 3 || isequal(SPM.Vbeta(1).dim, size(mask)), 'Different dimensions between mask and betas');

    % extract data and design matrix from confound GLM
    %

    Y = spm_data_read(SPM.xY.VY, find(mask)); % BOLD data

    X = SPM.xX.X; % original design matrix
    K = SPM.xX.K; % high-pass filter
    W = SPM.xX.W; % whitening matrix
    KWX = SPM.xX.xKXs.X; % high-pass filtered & whitened design matrix

    R = spm_sp('r',SPM.xX.xKXs); % residual forming matrix R = I - X * pinv(X)

    KWY = spm_filter(K,W*Y); % high-pass filtered & whitened data

    % convert filter K to matrix form
    % see spm_filter.m
    for s = 1:length(K)
        I(K(s).row,K(s).row) = eye(length(K(s).row));
        X0X0(K(s).row, K(s).row) = K(s).X0*K(s).X0';
    end
    K = I - X0X0; % high-pass filter matrix

    assert(immse(K*W*X, KWX) < 1e-15);
    assert(immse(K*W*Y, KWY) < 1e-15);
    assert(immse(R, eye(size(X,1)) - SPM.xX.xKXs.u*SPM.xX.xKXs.u') < 1e-15);
end
