%function fit_gp(subj, use_smooth, glmodel, mask, what, debug)

    subj = 1;
    use_smooth = true;
    glmodel = 21;
    %mask = 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii';
    mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
    what = 'theory';
    debug = true;

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end

    if ~exist('debug', 'var')
        debug = false;
    end

    assert(ismember(what, {'theory', 'sprite', 'interaction', 'termination'}));


    [~,maskname,~] = fileparts(mask);
    filename = sprintf('mat/fit_gp_HRR_subj=%d_us=%d_glm=%d_mask=%s_%s.mat', subj, use_smooth, glmodel, maskname, what);
    filename

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);

    addpath(genpath('/ncf/gershman/Lab/scripts/gpml'));


    fprintf('loading BOLD for subj %d\n', subj);
    tic
    [Y, K, W, R] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
    toc

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
    sigmas = logspace(-3, 4, 20);

    %{
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

        disp('    ... for CV');
        tic

        %
        % CV TODO dedupe
        %
        for p = 1:n_partitions % for each partition
            train = partition_id ~= p;
            n_CV = sum(train);

            % from infLikGauss.m
            sn2_CV(p,j) = exp(2*hyp.lik); W_CV = ones(n_CV,1)/sn2_CV(p,j);            % noise variance of likGauss
            K_gp_CV{p,j} = apx(hyp,covfun,x(train),[]);                        % set up covariance approximation
            [ldB2_CV(p,j),solveKiW_CV{p,j},~,~,~] = K_gp_CV{p,j}.fun(W_CV); % obtain functionality depending on W

            % compute inverse the standard way
            I_CV = eye(n_CV);
            invKi_CV{p,j} = (ker(train, train) + s^2 * I_CV)^(-1);
        end

        toc
    end

    %}

    fprintf('solving GP for subj %d, %d voxels\n', subj, size(Y,2));
    tic

    sigma = nan(1, size(Y,2)); % fitted noise std's
    loglik = nan(1, size(Y,2)); % marginal log likelihoods
    R2 = nan(1, size(Y,2)); % R^2
    adjR2 = nan(1, size(Y,2)); % adjusted R^2
    r = nan(1, size(Y,2)); % Pearson correlation

    % same but cross-validated
    sigma_CV = nan(n_partitions, size(Y,2));
    loglik_CV = nan(n_partitions, size(Y,2)); % PREDICTIVE log likelihood!
    R2_CV = nan(n_partitions, size(Y,2));
    adjR2_CV = nan(n_partitions, size(Y,2));
    r_CV = nan(n_partitions, size(Y,2));

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

        % for each sigma, compute marginal likelihood = loglik = -NLZ
        %
        for j = 1:length(sigmas)
            s = sigmas(j);
                
            % from infGaussLik
            nlz(j) = y'*invKi{j}*y/2 + ldB2(j) + n*log(2*pi*sn2(j))/2;    % -log marginal likelihood TODO there is no sn2 term in Eq 2.30 in the Rasmussen book?

            if debug
                % sanity checks
                hyp.lik = log(s);

                % GP log lik
                nlz_gp(j) = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);

                % from infGaussLik
                alpha = solveKiW{j}(y);
                nlz_gp2 = y'*alpha/2 + ldB2(j) + n*log(2*pi*sn2(j))/2;    % -log marginal likelihood

                assert(immse(nlz_gp(j), nlz_gp2) < 1e-3);
                assert(immse(nlz_gp(j), nlz(j)) < 1e-3);
            end
        end

        % pick best sigma
        [~,j] = min(nlz);
        sigma(i) = sigmas(j);
        loglik(i) = -nlz(j); % marginal log lik

        % posterior predictive on observed dataset (Eq. 2.23 in Rasmussen)
        % btw v slow... much slower than nlz computation
        y_hat = ker * invKi{j} * y;

        [R2(i), adjR2(i)] = calc_R2(y, y_hat, 1);
        r(i) = corr(y_hat, y);

        if debug
            hyp.lik = log(sigma(i));
            y_hat_gp = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y, x);

            assert(immse(y_hat, y_hat_gp) < 1e-10);

            %figure;
            %hold on;
            %plot(y);
            %plot(y_hat);
            %plot(y_hat_gp);
            %legend({'y', 'y_hat', 'y_hat_gp'}, 'interpreter', 'none');
        end

        %
        % same thing but cross-validated TODO dedupe
        %

        y_hat_CV = nan(size(y));
        for p = 1:n_partitions % for each partition
            train = partition_id ~= p;
            test = partition_id == p;

            n_CV = sum(train);

            % for each sigma, compute marginal likelihood = loglik = -NLZ
            %
            for j = 1:length(sigmas)
                s = sigmas(j);
                    
                % from infGaussLik
                nlz_CV(p,j) = y(train)'*invKi_CV{p,j}*y(train)/2 + ldB2_CV(p,j) + n_CV*log(2*pi*sn2_CV(p,j))/2;    % -log marginal likelihood TODO there is no sn2 term in Eq 2.30 in the Rasmussen book?

                if debug
                    % sanity checks
                    hyp.lik = log(s);

                    % GP log lik
                    nlz_gp_CV(p,j) = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train));

                    % from infGaussLik
                    alpha_CV = solveKiW_CV{p,j}(y(train));
                    nlz_gp2_CV = y(train)'*alpha_CV/2 + ldB2_CV(p,j) + n_CV*log(2*pi*sn2_CV(p,j))/2;    % -log marginal likelihood

                    assert(immse(nlz_gp_CV(p,j), nlz_gp2_CV) < 1e-3);
                    assert(immse(nlz_gp_CV(p,j), nlz_CV(p,j)) < 1e-3);
                end
            end

            % pick best sigma
            [~,j] = min(nlz_CV(p,:));
            sigma_CV(p,i) = sigmas(j);

            % compute predictive log lik (notice difference -- before it was marginal!)
            % TODO why are they all identical?
            [~,~,~,~,lp] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train), x(test), y(test));
            loglik_CV(p,i) = sum(lp);

            % posterior predictive on observed dataset (Eq. 2.23 in Rasmussen)
            y_hat_CV(test) = ker(test, train) * invKi_CV{p,j} * y(train);

            [R2_CV(p,i), adjR2_CV(p,i)] = calc_R2(y(test), y_hat(test), 1);
            r_CV(p,i) = corr(y_hat(test), y(test));

            if debug
                hyp.lik = log(sigma_CV(p,i));
                y_hat_gp_CV = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train), x(test));

                assert(immse(y_hat_CV(test), y_hat_gp_CV) < 1e-10);
            end
        end

        %figure;
        %hold on;
        %plot(y);
        %plot(y_hat_CV);
        %legend({'y', 'y_hat_CV'}, 'interpreter', 'none');

    end

    toc

    %{
    filename
    clear invKi
    clear invKi_gp
    clear I
    clear K
    clear K_gp
    clear L
    clear R
    clear hyp
    clear covhyp
    clear ker
    clear post
    clear Y
    save(filename, '-v7.3');
    %}

    disp('Done');

%end




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
