% meat of fit_gp_CV.m, helper for decode_gp_CV.m TODO dedupe w/ fit_gp_CV.m

function r_CV = fit_gp_CV_helper(subj, use_smooth, glmodel, mask, ker, sigma, fast, debug)

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end

    if ~exist('fast', 'var')
        fast = true; % no ceiling, sigma = 1 only, no logpredlik
    end
    if ~exist('debug', 'var')
        debug = false;
    end

    fprintf('loading BOLD for subj %d\n', subj);
    tic
    [Y, K, W, R, run_id] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
    toc

    assert(size(Y,1) == size(ker,1));
    assert(size(Y,1) == size(ker,2));

    % whiten, filter & project out nuisance regressors
    Y = R*K*W*Y;
    ker = R*K*W*ker*W'*K'*R';

    % get partitions from RSA 3
    rsa = vgdl_create_rsa(3, subj);
    partition_id = rsa.model(1).partitions;
    assert(size(partition_id, 1) == size(Y, 1));
    n_partitions = max(partition_id);

    % find nearest symmetric positive definite matrix (it's not b/c of numerical issues, floating points, etc.)
    ker = nearestSPD(ker);

    % init GP stuff
    %
    n = size(ker, 1);
    x = [1:n]';
    ceil_x = x + run_id * 10000; % do not extrapolate RBF across runs

    % init functions and  hyperparams from kernel
    %meanfun = {@meanDiscrete, n};
    meanfun = @meanConst;
    covfun = {@covDiscrete, n};
    likfun = @likGauss;
    hyp = hyp_from_ker(ker);


