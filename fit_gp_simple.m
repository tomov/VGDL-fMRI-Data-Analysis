% meat of fit_gp_CV.m, helper for decode_gp_CV.m TODO dedupe w/ fit_gp_CV.m and fit_gp_CV_simple.m
% unlike fit_gp_CV_simple, this one does'nt do the CV

function [nlz] = fit_gp_CV(Y, ker, x, y, meanfun, covfun, likfun, debug)

    %{
    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end
    %}

    if ~exist('debug', 'var')
        debug = false;
    end

    fast = true;

    addpath(genpath('/ncf/gershman/Lab/scripts/gpml'));

    % load mask
    %{
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);

    %fprintf('loading BOLD for subj %d\n', subj);
    [Y, K, W, R, run_id] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);

    assert(size(Y,1) == size(ker,1));
    assert(size(Y,1) == size(ker,2));

    % whiten, filter & project out nuisance regressors
    Y = R*K*W*Y;
    ker = R*K*W*ker*W'*K'*R';
    %}

    % find nearest symmetric positive definite matrix (it's not b/c of numerical issues, floating points, etc.)
    ker = nearestSPD(ker);
    hyp = hyp_from_ker(ker);

    % init GP stuff
    %
    %{
    n = size(ker, 1);
    x = [1:n]';
    ceil_x = x + run_id * 10000; % do not extrapolate RBF across runs

    % init functions and  hyperparams from kernel
    %meanfun = {@meanDiscrete, n};
    meanfun = @meanConst;
    covfun = {@covDiscrete, n};
    likfun = @likGauss;
    %}


    sigmas = 1; % TODO param


    % precompute (K + sigma^2 I) ^ (-1) for every sigma
    % also K(X*,X) * (K(X,X) + sigma^2 I) ^ (-1)  (predictive Eq. 2.23)
    %
    for j = 1:length(sigmas)
        s = sigmas(j);

        %fprintf('inverting kernel for subj %d, sigma %.3f\n', subj, s);

        [sn2(j), ...
         ldB2(j), ...
         solveKiW{j}, ...
         invKi{j}] = calc_invKi(ker, x, s, hyp, covfun);

        ker_invKi{j} = ker * invKi{j};

    end


    %assert(size(Y,2) == 1);

    %fprintf('solving GP for subj %d, %d voxels\n', subj, size(Y,2));

    sigma = nan(1, size(Y,2)); % fitted noise std's
    logmarglik = nan(1, size(Y,2)); % marginal log likelihoods
    logpredlik = nan(1, size(Y,2)); % predictive log likelihoods
    R2 = nan(1, size(Y,2)); % R^2
    adjR2 = nan(1, size(Y,2)); % adjusted R^2
    r = nan(1, size(Y,2)); % Pearson correlation
    MSE = nan(1, size(Y,2)); % MSE
    SMSE = nan(1, size(Y,2)); % SMSE


    % GP
    % for each voxel
    %
    for i = 1:size(Y, 2)

        y = Y(:,i);

        % no CV; use marginal likelihood for model comparison
        %
        train = logical(ones(size(y,1),1)); % all trials
        test = train; % test data = train data here

        [sigma(i), ...
         logmarglik(i), ...
         logpredlik(i), ...
         y_hat, ...
         R2(i), ...
         adjR2(i), ...
         r(i), ...
         MSE(i), ...
         SMSE(i)] = fit_gp_helper(x, y, train, test, ker, hyp, meanfun, covfun, likfun, sigmas, invKi, ldB2, sn2, solveKiW, ker_invKi, fast, debug);


        %{
        figure;
        hold on;
        plot(y);
        plot(y_hat);
        legend({'y', 'y_hat'}, 'interpreter', 'none');
        %}

    end


    nlz = -logmarglik;
