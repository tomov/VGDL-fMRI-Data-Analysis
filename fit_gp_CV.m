function fit_gp_CV(subj, use_smooth, glmodel, mask, model_name, what, fast, debug)

%{
    subj = 1;
    use_smooth = true;
    glmodel = 21;
    %mask = 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii';
    mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
    %mask = 'masks/mask_batchsize=1000_batch=2.nii';
    %mask = 'masks/mask.nii';
    model
    what = 'theory';
    fast = true;
    debug = false;
    %}

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



    [~,maskname,~] = fileparts(mask);
    if fast
        filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=%s_model=%s_%s_nsamples=1_fast.mat', subj, use_smooth, glmodel, maskname, model_name, what);
    else
        filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=%s_model=%s_%s_nsamples=1_notfast.mat', subj, use_smooth, glmodel, maskname, model_name, what);
    end
    filename

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);

    addpath(genpath('/ncf/gershman/Lab/scripts/gpml'));


    fprintf('loading BOLD for subj %d\n', subj);
    tic
    [Y, K, W, R, run_id] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
    toc

    fprintf('Memory usage: %.3f MB\n', monitor_memory_whos);

    fprintf('loading kernel for subj %d\n', subj);
    tic
    switch model_name
        case 'EMPA'
            assert(ismember(what, {'theory', 'sprite', 'interaction', 'termination'}));
            ker = load_HRR_kernel(subj, what);
        case 'DQN'
            assert(ismember(what, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2'}));
            ker = load_DQN_kernel(subj, what)
        case 'game'
            ker = load_game_kernel(EXPT, subj); % GLM 1 game id features
        otherwise
            assert(false, 'invalid model name')
    end
    toc

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
    ceil_x = x + run_id * 10000; % do not extrapolate RBF across runs

    % init functions and  hyperparams from kernel
    %meanfun = {@meanDiscrete, n};
    meanfun = @meanConst;
    covfun = {@covDiscrete, n};
    likfun = @likGauss;
    hyp = hyp_from_ker(ker);

    % same for noise ceiling model
    % from http://www.gaussianprocess.org/gpml/code/matlab/doc/, section 4a)
    ceil_meanfun = @meanConst;
    ceil_covfun = @covSEiso; % RBF
    ceil_likfun = @likGauss;
    ceil_hyp = struct('mean', 0, 'cov', [0; 0], 'lik', log(0.1));

    % grid search sigmas
    if fast
        sigmas = 1;
    else
        sigmas = logspace(-3, 3, 20);
    end

    % precompute (K + sigma^2 I) ^ (-1) for every sigma
    % also K(X*,X) * (K(X,X) + sigma^2 I) ^ (-1)  (predictive Eq. 2.23)
    %
    for j = 1:length(sigmas)
        s = sigmas(j);

        fprintf('inverting kernel for subj %d, sigma %.3f\n', subj, s);
        tic

        [sn2(j), ...
         ldB2(j), ...
         solveKiW{j}, ...
         invKi{j}] = calc_invKi(ker, x, s, hyp, covfun);

        ker_invKi{j} = ker * invKi{j};

        toc

        disp('    ... for CV');
        tic

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

        toc
    end




    fprintf('solving GP for subj %d, %d voxels\n', subj, size(Y,2));
    tic

    sigma = nan(1, size(Y,2)); % fitted noise std's
    logmarglik = nan(1, size(Y,2)); % marginal log likelihoods
    logpredlik = nan(1, size(Y,2)); % predictive log likelihoods
    R2 = nan(1, size(Y,2)); % R^2
    adjR2 = nan(1, size(Y,2)); % adjusted R^2
    r = nan(1, size(Y,2)); % Pearson correlation
    MSE = nan(1, size(Y,2)); % MSE
    SMSE = nan(1, size(Y,2)); % SMSE

    % same but cross-validated
    sigma_CV = nan(n_partitions, size(Y,2));
    logmarglik_CV = nan(n_partitions, size(Y,2));
    logpredlik_CV = nan(n_partitions, size(Y,2));
    R2_CV = nan(n_partitions, size(Y,2));
    adjR2_CV = nan(n_partitions, size(Y,2));
    r_CV = nan(n_partitions, size(Y,2));
    MSE_CV = nan(n_partitions, size(Y,2));
    SMSE_CV = nan(n_partitions, size(Y,2));

    % same but noise ceiling kernel
    ceil_sigma = nan(1, size(Y,2)); % fitted noise std's
    ceil_logmarglik = nan(1, size(Y,2)); % marginal log likelihoods
    ceil_logpseudolik = nan(1, size(Y,2)); % predictive log likelihoods
    ceil_R2 = nan(1, size(Y,2)); % R^2
    ceil_adjR2 = nan(1, size(Y,2)); % adjusted R^2
    ceil_r = nan(1, size(Y,2)); % Pearson correlation
    ceil_MSE = nan(1, size(Y,2)); % MSE
    ceil_SMSE = nan(1, size(Y,2)); % SMSE 

    % same but cross-validated
    ceil_sigma_CV = nan(n_partitions, size(Y,2));
    ceil_logmarglik_CV = nan(n_partitions, size(Y,2));
    ceil_logpseudolik_CV = nan(n_partitions, size(Y,2));
    ceil_R2_CV = nan(n_partitions, size(Y,2));
    ceil_adjR2_CV = nan(n_partitions, size(Y,2));
    ceil_r_CV = nan(n_partitions, size(Y,2));
    ceil_MSE_CV = nan(n_partitions, size(Y,2));
    ceil_SMSE_CV = nan(n_partitions, size(Y,2));

    % GP
    % for each voxel
    %
    for i = 1:size(Y, 2)

        if (fast && mod(i,1000) == 0) || (~fast)
            i
            toc
            tic
        end

        y = Y(:,i);

        % no CV; use marginal likelihood for model comparison
        %
        train = partition_id > 0; % all trials
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

        % fit noise ceiling model
        if ~fast
            [ceil_sigma(i), ...
             ceil_logmarglik(i), ...
             ceil_logpseudolik(i), ...
             ceil_y_hat, ...
             ceil_R2(i), ...
             ceil_adjR2(i), ...
             ceil_r(i), ...
             ceil_MSE(i), ...
             ceil_SMSE(i), ...
             ceil_hyp] = minimize_gp_helper(ceil_x, y, train, test, ceil_hyp, ceil_meanfun, ceil_covfun, ceil_likfun, debug);
        end

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

            % noise ceiling model
            if ~fast
                [ceil_sigma_CV(p,i), ...
                 ceil_logmarglik_CV(p,i), ...
                 ceil_logpseudolik_CV(p,i), ...
                 ceil_y_hat_CV(test), ...
                 ceil_R2_CV(p,i), ...
                 ceil_adjR2_CV(p,i), ...
                 ceil_r_CV(p,i), ...
                 ceil_MSE_CV(p,i), ...
                 ceil_SMSE_CV(p,i), ...
                 ceil_hyp_CV] = minimize_gp_helper(ceil_x, y, train, test, ceil_hyp, ceil_meanfun, ceil_covfun, ceil_likfun, debug);
            end

        end

        %figure;
        %hold on;
        %plot(y);
        %plot(y_hat);
        %legend({'y', 'y_hat'}, 'interpreter', 'none');

        %figure;
        %hold on;
        %plot(y);
        %plot(ceil_y_hat_CV);
        %legend({'y', 'ceil_y_hat', 'ceil_y_hat_CV'}, 'interpreter', 'none');

    end

    toc

    fprintf('Memory usage: %.3f MB\n', monitor_memory_whos);

    filename
    %{
    clear invKi_CV
    clear null_invKi_CV
    clear invKi
    clear null_invKi
    clear invKi_gp
    clear I
    clear K
    clear K_gp
    clear L
    clear R
    clear hyp
    clear null_hyp
    clear covhyp
    clear ker
    clear null_ker
    clear post
    clear Y
    %}
    save(filename, 'sigma', 'logmarglik', 'logpredlik', 'R2', 'adjR2', 'r', 'MSE', 'SMSE', ... 
                   'sigma_CV', 'logmarglik_CV', 'logpredlik_CV', 'R2_CV', 'adjR2_CV', 'r_CV', 'MSE_CV', 'SMSE_CV', ...
                   'ceil_sigma', 'ceil_logmarglik', 'ceil_logpseudolik', 'ceil_R2', 'ceil_adjR2', 'ceil_r', 'ceil_MSE', 'ceil_SMSE', ...
                   'ceil_sigma_CV', 'ceil_logmarglik_CV', 'ceil_logpseudolik_CV', 'ceil_R2_CV', 'ceil_adjR2_CV', 'ceil_r_CV', 'ceil_MSE_CV', 'ceil_SMSE_CV', ...
                   'subj', 'use_smooth', 'glmodel', 'mask', 'what', 'sigmas', 'n', ...
    '-v7.3');

    disp('Done');

end





% fit gp using minimize(), for noise ceiling model
% use LOO-CV to get predictions on test data
%
function [sigma, logmarglik, logpseudolik, y_hat, R2, adjR2, r, MSE, SMSE, hyp] = minimize_gp_helper(x, y, train, test, hyp, meanfun, covfun, likfun, debug)

    % fit mean, kernel, and sigma; 25 iters
    hyp = minimize(hyp, @gp, -25, @infGaussLik, meanfun, covfun, likfun, x(train), y(train));
    sigma = exp(hyp.lik);

    % compute marginal lik on training data
    nlz = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train));
    logmarglik = -nlz;

    % compute LOO likelihood = pseudolikelihood on test data (Eq. 5.12 in Rasmussen)
    % call it manually to get predictions
    [~, nloo, ~, loo] = infLOO(hyp, {meanfun}, {covfun}, {likfun}, x(test), y(test));
    logpseudolik = -nloo;
    y_hat = loo.ymu;

    if debug
        % call infLOO via gp()
        nloo2 = gp(hyp, @infLOO, meanfun, covfun, likfun, x(test), y(test));
        assert(immse(nloo2, nloo) < 1e-20);
        assert(immse(nloo, -sum(loo.lp)) < 1e-20);
        assert(immse(loo.ymu, loo.fmu) < 1e-20);

        % actually compute LOO manually (SUPER SLOW)
        ix = find(test);
        y_hat_2 = nan(size(y(test)));
        lp = nan(size(y(test)));
        for i = 1:length(ix)
            test(ix(i)) = 0;
            [y_hat_2(i), ~, ~, ~, lp(i)] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(test), y(test), x(ix(i)), y(ix(i)));
            test(ix(i)) = 1;
        end

        assert(immse(y_hat_2, loo.ymu) < 1e-10);
        assert(immse(lp, loo.lp) < 1e-10);
    end

    % R^2
    assert(all(train == test) || ~any(train == test));
    if all(train == test)
        % # of free (hyper)parameters
        p = numel(hyp.cov) + numel(hyp.lik);
        if isfield(hyp, 'mean')
            p = p + numel(hyp.mean);
        end
    else
        p = 0; % 0 params, b/c evaluating on held out data
    end
    [R2, adjR2] = calc_R2(y(test), y_hat, p);

    % Pearson
    r = corr(y_hat, y(test));

    % MSE and SMSE, Sec 2.5 in Rasmussen
    MSE = immse(y(test), y_hat);
    SMSE = MSE / var(y(test));
end


% load game identity kernel based on GLM 1
%
function [ker] = load_game_kernel(EXPT, subj_id)
    feature_glm = 1;

    modeldir = fullfile(EXPT.modeldir,['model',num2str(feature_glm)],['subj',num2str(subj_id)]);
    load(fullfile(modeldir,'SPM.mat'));

    if subj_id <= 11
        game_names = {'vgfmri3_chase','vgfmri3_helper','vgfmri3_bait','vgfmri3_lemmings','vgfmri3_plaqueAttack','vgfmri3_zelda'};
    else:
        game_names = {'vgfmri4_chase', 'vgfmri4_helper', 'vgfmri4_bait', 'vgfmri4_lemmings', 'vgfmri4_avoidgeorge', 'vgfmri4_zelda'};
    end

    % use convolved game boxcars from GLM 1
    % that will take HRF into account and even the slight overlap between games
    features = [];
    for g = 1:length(game_names)
        which = contains(SPM.xX.name, game_names(g));
        assert(sum(which) == 3);
        feature = sum(SPM.xX.X(:,which), 2); % merge game boxcars from different runs
        features(:,g) = feature;
    end

    % see gen_subject_kernels() in HRR.py

    sigma_w = 1; % TODO fit; matching with sigma_w in HRR.py

    Sigma_w = eye(size(features,2)) * sigma_w; % Sigma_p in Rasmussen, Eq. 2.4

    ker = features * Sigma_w * features'; % K in Rasmussen, Eq. 2.12
end


% load HRR kernel
%
function [ker] = load_HRR_kernel(subj_id, what)
    filename = sprintf('mat/HRR_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=1_sigma_w=1.000_norm=1.mat', subj_id);
    load(filename, 'theory_kernel', 'sprite_kernel', 'interaction_kernel', 'termination_kernel');

    ker = eval([what, '_kernel']);
end


function [ker] = load_DQN_kernel(subj_id, what)
    filename = sprintf('mat/DQN_subject_kernel_subj=%d_sigma_w=1.000_norm=1.mat', subj_id);
    load(filename, 'layer_conv1_output_kernel', 'layer_conv2_output_kernel', 'layer_conv3_output_kernel', 'layer_linear1_output_kernel', 'layer_linear2_output_kernel');

    ker = eval(['layer_', what, '_output_kernel']);
end
