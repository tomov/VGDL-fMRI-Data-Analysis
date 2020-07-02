clear all;
close all;

EXPT = vgdl_expt();
glmodel = 21; % control GLM
mask = 'masks/ROI_x=48_y=12_z=30_85voxels_Sphere6.nii';
subjects = 1:length(EXPT.subject);

debug = true;


% load mask
[mask_format, mask, Vmask] = get_mask_format_helper(mask);

addpath(genpath('/ncf/gershman/Lab/scripts/gpml'));


subj = subjects(1); % TODO for each subject

fprintf('loading BOLD for subj %d\n', subj);
tic
[Y, K, W, R] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
toc

fprintf('loading kernel for subj %d\n', subj);
tic
ker = load_kernel(subj);
toc

% whiten, filter & project out nuisance regressors
Y = R*K*W*Y;
ker = R*K*W*ker*W'*K'*R';



%{
rng(334);

sigma = 0.01;
tau = 3;
X = [ones(50,1) rand(50,3)];
b = [0; normrnd(0, tau, 3, 1)];
e = normrnd(0, sigma, size(X,1), 1);
y = X * b + e;
Sigma_p = tau.^2 * eye(size(X,2));



Y = [y y y];
ker = X * Sigma_p * X';
%}




%
% GP shenanighans
%
% naively, we would just loop over voxels, and run minimize() with gp() on each voxel, fit the sigma (noise std), and get the LME 
% however, we have two problems:
% 1) this will compute the inverse of (K + sigma^2 I) separately for every voxel (and every sigma) thus duplicating a lot of compute
% 2) to ensure only sigma is fit (and not the whole kernel and mean function), we have to impose hyperpriors (see http://www.gaussianprocess.org/gpml/code/matlab/doc/), however that is suuuper slow (even though we are optimizing way fewer parameters!)
%
% Therefore, we precompute the inverse of (K + sigma^2 I) for a grid of sigmas, then reuse it across voxels (note that this might be numerically unstable -- see A4 in the Rasmussen GP book).
% 1) is fixed because we now only compute the inverse a fixed # of times, not for each voxel
% 2) is fixed beacuse we're doing a grid search manually, rather than relying on the fancy minimize() (we don't need the exact sigmas anyway, ballparks are ok)
%
% note that their computation is different from the one in test_gp b/c 1) the inverse is approximate, and 2) the log determinant is approximate
%

% init GP stuff
%
n = size(ker, 1);
x = [1:n]';

meanfun = {@meanDiscrete, n};
covfun = {@covDiscrete, n};
likfun = @likGauss;

ker = nearestSPD(ker);  % find nearest symmetric positive definite matrix (it's not b/c of numerical issues, floating points, etc.)
L = cholcov(ker); % not chol b/c sometimes positive semi-definite
L(1:(n+1):end) = log(diag(L));  % see covDiscrete.m
covhyp = L(triu(true(n)));  % see covDiscrete.m

% GP hyperparams
hyp = struct('mean', zeros(1, n), 'cov', covhyp, 'lik', nan);

% grid search sigmas
sigmas = logspace(-10, 10, 21);

% precompute (K + sigma^2 I) ^ (-1) for every sigma
%
for j = 1:length(sigmas)
    s = sigmas(j);
    hyp.lik = log(s);

    fprintf('inverting kernel for subj %d, sigma %.3f\n', subj, s);
    tic

    % from infLikGauss.m
    sn2(j) = exp(2*hyp.lik); W = ones(n,1)/sn2(j);            % noise variance of likGauss
    K_gp{j} = apx(hyp,covfun,x,[]);                        % set up covariance approximation
    [ldB2(j),solveKiW{j},dW,dhyp,post.L] = K_gp{j}.fun(W); % obtain functionality depending on W

    assert(immse(K_gp{j}.mvm(eye(size(ker))), ker) < 1e-15); % should be identical

    % compute inverse the gp()-way (for sanity check)
    invKi_gp = solveKiW{j}(eye(n));

    % compute inverse the standard way
    I = eye(n);
    invKi{j} = (ker + s^2 * I)^(-1);

    %assert(immse(invKi{j}, invKi_gp) < 1e-15); % those are different! what matters is if the log likelihoods below match

    toc
end


% GP
% for each voxel
%
for i = 1:size(Y, 2)
    fprintf('solving GP for subj %d, voxel %d\n', subj, i);
    tic

    y = Y(:,i);

    % for each sigma, compute marginal likelihood = loglik = -NLZ
    %
    for j = 1:length(sigmas)
        s = sigmas(j);
            
        % from infGaussLik
        nlz(j) = y'*invKi{j}*y + ldB2(j) + n*log(2*pi*sn2(j))/2;    % -log marginal likelihood TODO there is no sn2 term in Eq 2.30 in the Rasmussen book?

        if debug
            % sanity checks
            hyp.lik = log(s);

            % GP log lik
            nlz_gp(j) = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);

            % from infGaussLik
            alpha = solveKiW{j}(y);
            nlz_gp2 = y'*alpha/2 + ldB2(j) + n*log(2*pi*sn2(j))/2;    % -log marginal likelihood

            assert(immse(nlz_gp(j), nlz_gp2) < 1e-15);
            %assert(immse(nlz_gp(j), nlz(j)) < 1e-15); % not equal.... bummer
        end
    end

    [~,j] = min(nlz);
    sigma(i) = sigmas(j);
    loglik(i) = -nlz(j);

    if debug
        % make sure at least sigmas are the same
        [~,j1] = min(nlz_gp);
        assert(j == j1);
    end

    toc
end




function [ker] = load_kernel(subj_id)
    filename = sprintf('mat/HRR_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=10_sigma_w=1.000_norm=1.mat', subj_id);
    load(filename, 'theory_kernel', 'sprite_kernel', 'interaction_kernel', 'termination_kernel');

    ker = theory_kernel;
end


function [Y, K, W, R] = load_BOLD(EXPT, glmodel, subj_id, mask, Vmask)
    % load subject data
    % Y = raw BOLD
    % K = filter matrix
    % W = whitening matrix
    % R = residual forming matrix
    %
    % 
    % Y = Xb + Gg + e              (Y = BOLD, X = HRR features; G = nuisance regressors (X in GLM 21), e = noise)
    % KWY = KWXb + KWGg + KWe      (whiten & filter w.r.t. G, as in GLM 21)
    % R = (I - KWG (KWG)^+)        (residual forming matrix w.r.t. G from GLM 21)
    % RKWY = RKWXb + 0 + RKWe      (project out nuisance regressors)
    % but...use GP regression instead
    %
    % TODO RKWe not gaussian, also maybe correlated => need to whiten again after GP fit

    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj_id)]);
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
