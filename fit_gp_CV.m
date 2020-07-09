function fit_gp_CV(subj, use_smooth, glmodel, mask, what, debug)

%{
    subj = 1;
    use_smooth = true;
    glmodel = 21;
    %mask = 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii';
    mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
    %mask = 'masks/mask_batchsize=1000_batch=2.nii';
    what = 'theory';
    debug = true;
    %}

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end

    if ~exist('debug', 'var')
        debug = true;
    end

    assert(ismember(what, {'theory', 'sprite', 'interaction', 'termination'}));


    [~,maskname,~] = fileparts(mask);
    filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=%s_%s.mat', subj, use_smooth, glmodel, maskname, what);
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
    ker = load_HRR_kernel(subj, what);
    null_ker = load_null_kernel(EXPT, subj); % GLM 1 game id features
    toc

    % whiten, filter & project out nuisance regressors
    Y = R*K*W*Y;
    ker = R*K*W*ker*W'*K'*R';
    null_ker = R*K*W*null_ker*W'*K'*R'; % TODO technically, this should project out the null kernel! b/c game identity features were already part of GLM 1

    % get partitions from RSA 3
    rsa = vgdl_create_rsa(3, subj);
    partition_id = rsa.model(1).partitions;
    assert(size(partition_id, 1) == size(Y, 1));
    n_partitions = max(partition_id);

    % find nearest symmetric positive definite matrix (it's not b/c of numerical issues, floating points, etc.)
    ker = nearestSPD(ker);
    null_ker = nearestSPD(null_ker);


    %{
    % for sanity checks
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
    ceil_x = x + run_id * 10000; % do not extrapolate RBF across runs

    % init functions and  hyperparams from kernel
    %meanfun = {@meanDiscrete, n};
    meanfun = @meanConst;
    covfun = {@covDiscrete, n};
    likfun = @likGauss;
    hyp = hyp_from_ker(ker);

    % same for null model
    null_meanfun = @meanConst;
    null_covfun = {@covDiscrete, n};
    null_likfun = @likGauss;
    null_hyp = hyp_from_ker(null_ker);

    % same for noise ceiling model
    % from http://www.gaussianprocess.org/gpml/code/matlab/doc/, section 4a)
    ceil_meanfun = [];
    ceil_covfun = @covSEiso; % RBF
    ceil_likfun = @likGauss;
    ceil_hyp = struct('cov', [0; 0], 'lik', log(0.1));

    % grid search sigmas
    %sigmas = logspace(-3, 4, 20);
    sigmas = 1;

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

        [null_sn2(j), ...
         null_ldB2(j), ...
         null_solveKiW{j}, ...
         null_invKi{j}] = calc_invKi(null_ker, x, s, null_hyp, covfun);

        null_ker_invKi{j} = null_ker * null_invKi{j};

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

            [null_sn2_CV(p,j), ...
             null_ldB2_CV(p,j), ...
             null_solveKiW_CV{p,j}, ...
             null_invKi_CV{p,j}] = calc_invKi(null_ker(train, train), x(train), s, null_hyp, covfun);

            null_ker_invKi_CV{p,j} = null_ker(test, train) * null_invKi_CV{p,j};
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

    % same but null kernel
    null_sigma = nan(1, size(Y,2)); % fitted noise std's
    null_logmarglik = nan(1, size(Y,2)); % marginal log likelihoods
    null_logpredlik = nan(1, size(Y,2)); % predictive log likelihoods
    null_R2 = nan(1, size(Y,2)); % R^2
    null_adjR2 = nan(1, size(Y,2)); % adjusted R^2
    null_r = nan(1, size(Y,2)); % Pearson correlation
    null_MSE = nan(1, size(Y,2)); % MSE 
    null_SMSE = nan(1, size(Y,2)); % SMSE 

    % same but cross-validated
    null_sigma_CV = nan(n_partitions, size(Y,2));
    null_logmarglik_CV = nan(n_partitions, size(Y,2));
    null_logpredlik_CV = nan(n_partitions, size(Y,2));
    null_R2_CV = nan(n_partitions, size(Y,2));
    null_adjR2_CV = nan(n_partitions, size(Y,2));
    null_r_CV = nan(n_partitions, size(Y,2));
    null_MSE_CV = nan(n_partitions, size(Y,2));
    null_SMSE_CV = nan(n_partitions, size(Y,2));

    % same but noise ceiling kernel
    ceil_sigma = nan(1, size(Y,2)); % fitted noise std's
    ceil_logmarglik = nan(1, size(Y,2)); % marginal log likelihoods
    ceil_logpredlik = nan(1, size(Y,2)); % predictive log likelihoods
    ceil_R2 = nan(1, size(Y,2)); % R^2
    ceil_adjR2 = nan(1, size(Y,2)); % adjusted R^2
    ceil_r = nan(1, size(Y,2)); % Pearson correlation
    ceil_MSE = nan(1, size(Y,2)); % MSE
    ceil_SMSE = nan(1, size(Y,2)); % SMSE 

    % same but cross-validated
    ceil_sigma_CV = nan(n_partitions, size(Y,2));
    ceil_logmarglik_CV = nan(n_partitions, size(Y,2));
    ceil_logpredlik_CV = nan(n_partitions, size(Y,2));
    ceil_R2_CV = nan(n_partitions, size(Y,2));
    ceil_adjR2_CV = nan(n_partitions, size(Y,2));
    ceil_r_CV = nan(n_partitions, size(Y,2));
    ceil_MSE_CV = nan(n_partitions, size(Y,2));
    ceil_SMSE_CV = nan(n_partitions, size(Y,2));

    % GP
    % for each voxel
    %
    for i = 1:size(Y, 2)
        %fprintf('solving GP for subj %d, voxel %d\n', subj, i);
        %tic

        if mod(i,1000) == 0
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
         SMSE(i)] = fit_gp_helper(x, y, train, test, ker, hyp, meanfun, covfun, likfun, sigmas, invKi, ldB2, sn2, solveKiW, ker_invKi, debug);

        % fit null model
        [null_sigma(i), ...
         null_logmarglik(i), ...
         null_logpredlik(i), ...
         null_y_hat, ...
         null_R2(i), ...
         null_adjR2(i), ...
         null_r(i), ...
         null_MSE(i), ... 
         null_SMSE(i)] = fit_gp_helper(x, y, train, test, null_ker, null_hyp, null_meanfun, null_covfun, null_likfun, sigmas, null_invKi, null_ldB2, null_sn2, null_solveKiW, null_ker_invKi, debug);

         %{
        % fit noise ceiling model
        [ceil_sigma(i), ...
         ceil_logmarglik(i), ...
         ceil_logpredlik(i), ...
         ceil_y_hat, ...
         ceil_R2(i), ...
         ceil_adjR2(i), ...
         ceil_r(i), ...
         ceil_hyp] = minimize_gp_helper(ceil_x, y, train, test, ceil_hyp, ceil_meanfun, ceil_covfun, ceil_likfun);
         %}

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
             SMSE_CV(p,i)] = fit_gp_helper(x, y, train, test, ker, hyp, meanfun, covfun, likfun, sigmas, invKi_CV(p,:), ldB2_CV(p,:), sn2_CV(p,:), solveKiW_CV(p,:), ker_invKi_CV(p,:), debug);

            % null model
            [null_sigma_CV(p,i), ...
             null_logmarglik_CV(p,i), ...
             null_logpredlik_CV(p,i), ...
             null_y_hat_CV(test), ...
             null_R2_CV(p,i), ...
             null_adjR2_CV(p,i), ...
             null_r_CV(p,i), ...
             null_MSE_CV(p,i), ...
             null_SMSE_CV(p,i)] = fit_gp_helper(x, y, train, test, null_ker, null_hyp, null_meanfun, null_covfun, null_likfun, sigmas, null_invKi_CV(p,:), null_ldB2_CV(p,:), null_sn2_CV(p,:), null_solveKiW_CV(p,:), null_ker_invKi_CV(p,:), debug);

            % noise ceiling model
            %{
            [ceil_sigma_CV(p,i), ...
             ceil_logmarglik_CV(p,i), ...
             ceil_logpredlik_CV(p,i), ...
             ceil_y_hat_CV(test), ...
             ceil_R2_CV(p,i), ...
             ceil_adjR2_CV(p,i), ...
             ceil_r_CV(p,i), ...
             ceil_hyp_CV] = minimize_gp_helper(ceil_x, y, train, test, ceil_hyp, ceil_meanfun, ceil_covfun, ceil_likfun);
             %}

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
                   'null_sigma', 'null_logmarglik', 'null_logpredlik', 'null_R2', 'null_adjR2', 'null_r', 'null_MSE', 'null_SMSE', ...
                   'null_sigma_CV', 'null_logmarglik_CV', 'null_logpredlik_CV', 'null_R2_CV', 'null_adjR2_CV', 'null_r_CV', 'null_MSE_CV', 'null_SMSE_CV', ...
                   'ceil_sigma', 'ceil_logmarglik', 'ceil_logpredlik', 'ceil_R2', 'ceil_adjR2', 'ceil_r', 'ceil_MSE', 'ceil_SMSE', ...
                   'ceil_sigma_CV', 'ceil_logmarglik_CV', 'ceil_logpredlik_CV', 'ceil_R2_CV', 'ceil_adjR2_CV', 'ceil_r_CV', 'ceil_MSE_CV', 'ceil_SMSE_CV', ...
                   'subj', 'use_smooth', 'glmodel', 'mask', 'what', 'sigmas', ...
    '-v7.3');

    disp('Done');

end





% fit gp using minimize(), for noise ceiling model
%
function [sigma, logmarglik, logpredlik, y_hat, R2, adjR2, r, ceil_hyp, MSE, SMSE] = minimize_gp_helper(x, y, train, test, hyp, meanfun, covfun, likfun)

    % fit mean, kernel, and sigma; 100 iters
    ceil_hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfun, covfun, likfun, x(train), y(train));

    % compute marginal lik on training data
    nlz = gp(ceil_hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train));

    % compute predictive lik on observed test data
    % LOO-CV on test data for prediction
    % notice we're ignoring training data b/c it's from different runs
    % TODO too slow...
    ix = find(test);
    y_hat = nan(size(y(test)));
    lp = nan(size(y(test)));
    for i = 1:length(ix)
        test(ix(i)) = 0;
        [y_hat(i), ~, ~, ~, lp(i)] = gp(ceil_hyp, @infGaussLik, meanfun, covfun, likfun, x(test), y(test), x(ix(i)), y(ix(i)));
        test(ix(i)) = 1;
    end
    %[y_hat, ~, ~, ~, lp] = gp(ceil_hyp, @infGaussLik, meanfun, covfun, likfun, x(test), y(test), x(test), y(test));   -- don't do this; this double dips the test data and gives overly optimistic predictions

    sigma = exp(ceil_hyp.lik);
    logmarglik = -nlz;
    logpredlik = sum(lp);

    % R^2
    assert(all(train == test) || ~any(train == test));
    if all(train == test)
        % # of free (hyper)parameters
        p = numel(ceil_hyp.cov) + numel(ceil_hyp.lik);
        if isfield(ceil_hyp, 'mean')
            p = p + numel(ceil_hyp.mean);
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


% fit gp sigmas manually; optionally sanity check with gp() library
%
function [sigma, logmarglik, logpredlik, y_hat, R2, adjR2, r, MSE, SMSE] = fit_gp_helper(x, y, train, test, ker, hyp, meanfun, covfun, likfun, sigmas, invKi, ldB2, sn2, solveKiW, ker_invKi, debug)

    n = sum(train);

    % for each sigma, compute marginal likelihood = loglik = -NLZ
    %
    for j = 1:length(sigmas)
        s = sigmas(j);
            
        % from infGaussLik
        nlz(j) = y(train)'*invKi{j}*y(train)/2 + ldB2(j) + n*log(2*pi*sn2(j))/2;    % -log marginal likelihood TODO there is no sn2 term in Eq 2.30 in the Rasmussen book?

        if debug
            % sanity checks
            hyp.lik = log(s);

            % GP log lik
            nlz_gp(j) = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train));

            % from infGaussLik
            alpha = solveKiW{j}(y(train));
            nlz_gp2 = y(train)'*alpha/2 + ldB2(j) + n*log(2*pi*sn2(j))/2;    % -log marginal likelihood

            %{
            nlz_gp(j)
            nlz_gp2
            nlz(j)
            immse(nlz_gp(j), nlz_gp2)
            immse(nlz_gp(j), nlz(j))
            %}
            assert(immse(nlz_gp(j), nlz_gp2) < 1e-1);
            assert(immse(nlz_gp(j), nlz(j)) < 1e-1);
        end
    end

    % pick best sigma
    [~,j] = min(nlz);
    sigma = sigmas(j);
    logmarglik = -nlz(j); % marginal log lik

    % compute predictive log lik
    hyp.lik = log(sigma);
    %[~,~,~,~,lp] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train), x(test), y(test));
    %logpredlik = sum(lp);
    logpredlik = nan;

    % posterior predictive on test dataset (Eq. 2.23 in Rasmussen)
    y_hat = ker_invKi{j} * y(train);

    % R^2
    assert(all(train == test) || ~any(train == test));
    if all(train == test)
        p = 1; % # 1 (hyper)param = sigma; b/c we fit on same data 
    else
        p = 0; % 0 params, b/c evaluating on held out data
    end
    [R2, adjR2] = calc_R2(y(test), y_hat, 1);

    % Pearson
    r = corr(y_hat, y(test));

    % MSE and SMSE, Sec. 2.5 in Rasmussen
    MSE = immse(y(test), y_hat);
    SMSE = MSE / var(y(test));

    if debug
        hyp.lik = log(sigma);
        y_hat_gp = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train), x(test));

        assert(immse(y_hat, y_hat_gp) < 1e-10);

        y_hat_2 = ker(test, train) * invKi{j} * y(train); % too slow => precompute
        assert(immse(y_hat, y_hat_2) < 1e-10);


        %figure;
        %hold on;
        %plot(y);
        %plot(y_hat);
        %plot(y_hat_gp);
        %legend({'y', 'y_hat', 'y_hat_gp'}, 'interpreter', 'none');
    end

end


% load game identity kernel based on GLM 1
%
function [ker] = load_null_kernel(EXPT, subj_id)
    feature_glm = 1;

    modeldir = fullfile(EXPT.modeldir,['model',num2str(feature_glm)],['subj',num2str(subj_id)]);
    load(fullfile(modeldir,'SPM.mat'));

    game_names = {'vgfmri3_chase','vgfmri3_helper','vgfmri3_bait','vgfmri3_lemmings','vgfmri3_plaqueAttack','vgfmri3_zelda'};

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
    filename = sprintf('mat/HRR_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=10_sigma_w=1.000_norm=1.mat', subj_id);
    load(filename, 'theory_kernel', 'sprite_kernel', 'interaction_kernel', 'termination_kernel');

    ker = eval([what, '_kernel']);
end


function [Y, K, W, R, run_id] = load_BOLD(EXPT, glmodel, subj_id, mask, Vmask)
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

    for r = 1:length(SPM.Sess)
        run_id(SPM.Sess(r).row,:) = r;
    end
end



function [R2, adjR2] = calc_R2(y, y_hat, p)
    % calculate R^2 https://en.wikipedia.org/wiki/Coefficient_of_determination
    %
    % p = # of free parameters (e.g. sigma)
    %
    n = length(y);

    SStot = sum((y - mean(y)).^2); 
    assert(immse(SStot, var(y) * (n - 1)) < 1e-15);
    SSres = sum((y - y_hat).^2);

    R2 = 1 - SSres/SStot;

    adjR2 = 1 - (1 - R2) * (n - 1) / (n - p - 1);
end


% init hyperparams from kernel for gp() toolbox
%
function [hyp] = hyp_from_ker(ker)
    n = size(ker, 1);

    % init hyperparams for kernel
    L = cholcov(ker); % not chol b/c sometimes positive semi-definite
    L(1:(n+1):end) = log(diag(L));  % see covDiscrete.m
    covhyp = L(triu(true(n)));  % see covDiscrete.m

    % GP hyperparams
    %hyp = struct('mean', zeros(1, n), 'cov', covhyp, 'lik', nan);
    hyp = struct('mean', 0, 'cov', covhyp, 'lik', nan);
end


function [sn2, ldB2, solveKiW, invKi] = calc_invKi(ker, x, s, hyp, covfun)
    n = size(ker, 1);
    hyp.lik = log(s);

    % from infLikGauss.m
    sn2 = exp(2*hyp.lik); W = ones(n,1)/sn2;            % noise variance of likGauss
    K_gp = apx(hyp,covfun,x,[]);                        % set up covariance approximation
    [ldB2,solveKiW,dW,dhyp,post.L] = K_gp.fun(W); % obtain functionality depending on W

    assert(immse(K_gp.mvm(eye(size(ker))), ker) < 1e-15); % should be identical

    % compute inverse the gp()-way (for sanity check)
    invKi_gp = solveKiW(eye(n));

    % compute inverse the standard way
    I = eye(n);
    invKi = (ker + s^2 * I)^(-1);

    %assert(immse(invKi, invKi_gp) < 1e-15); % those are different! what matters is if the log likelihoods below match

end
