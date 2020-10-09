% meat of fit_gp_CV.m, helper for decode_gp_CV.m TODO dedupe w/ fit_gp_CV.m

%function [r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj, use_smooth, glmodel, mask, ker, debug)
function [r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj, use_smooth, glmodel, Y, ker, run_id, x, y, meanfun, covfun, likfun, partition_id, debug)

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end

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

    if ~exist('partition_id', 'var')
        % get partitions from RSA 3
        rsa = vgdl_create_rsa(3, subj);
        partition_id = rsa.model(1).partitions;
    end
    assert(size(partition_id, 1) == size(Y, 1));
    n_partitions = max(partition_id);

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

        %
        % CV
        %
        for p = 1:n_partitions % for each partition
            train = partition_id ~= p;
            test = partition_id == p;

            [sn2_CV(p,j), ...
             ldB2_CV(p,j), ...
             solveKiW_CV{p,j}, ...
             invKi_CV{p,j}] = calc_invKi(ker(train, train), x(train), s, hyp, covfun);

            ker_invKi_CV{p,j} = ker(test, train) * invKi_CV{p,j};
        end

    end


    %assert(size(Y,2) == 1);

    %fprintf('solving GP for subj %d, %d voxels\n', subj, size(Y,2));

    % same but cross-validated
    sigma_CV = nan(n_partitions, size(Y,2));
    logmarglik_CV = nan(n_partitions, size(Y,2));
    logpredlik_CV = nan(n_partitions, size(Y,2));
    R2_CV = nan(n_partitions, size(Y,2));
    adjR2_CV = nan(n_partitions, size(Y,2));
    r_CV = nan(n_partitions, size(Y,2));
    MSE_CV = nan(n_partitions, size(Y,2));
    SMSE_CV = nan(n_partitions, size(Y,2));


    % GP
    % for each voxel
    %
    for i = 1:size(Y, 2)

        y = Y(:,i);

        % CV; use predictive likelihood for model comparison
        %
        y_hat_CV = nan(size(y));
        for p = 1:n_partitions % for each partition
            train = partition_id ~= p;
            test = partition_id == p;

            [sigma_CV(p,i), ...
             logmarglik_CV(p,i), ...
             logpredlik_CV(p,i), ...
             y_hat_CV(test), ...
             R2_CV(p,i), ...
             adjR2_CV(p,i), ...
             r_CV(p,i), ...
             MSE_CV(p,i), ...
             SMSE_CV(p,i)] = fit_gp_helper(x, y, train, test, ker, hyp, meanfun, covfun, likfun, sigmas, invKi_CV(p,:), ldB2_CV(p,:), sn2_CV(p,:), solveKiW_CV(p,:), ker_invKi_CV(p,:), fast, debug);

        end

        %{
        figure;
        hold on;
        plot(y);
        plot(y_hat);
        legend({'y', 'y_hat'}, 'interpreter', 'none');
        %}

    end

