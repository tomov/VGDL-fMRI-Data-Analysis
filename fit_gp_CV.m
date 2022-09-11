function fit_gp_CV(subj, use_smooth, glmodel, mask, model_name, what, project, normalize, concat, novelty, fast, save_Y_hat, which_partitions, which_games, debug)

    %{
    clear all;
    subj = 1;
    use_smooth = true;
    glmodel = 1;
    %mask = 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii';
    %mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
    mask = 'masks/ROI_x=16_y=-94_z=22_1voxels_Sphere1.nii';
    %mask = 'masks/mask_batchsize=1000_batch=2.nii';
    %mask = 'masks/mask.nii';
    model_name = 'state'
    what = '';
    fast = true;
    project = true;
    debug = false;
    save_Y_hat = true;
    %}

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end

    if ~exist('fast', 'var')
        fast = true; % no ceiling, sigma = 1 only, no logpredlik
    end
    if ~exist('save_Y_hat', 'var')
        save_Y_hat = false;
    end
    if ~exist('debug', 'var')
        debug = false;
    end
    if ~exist('project', 'var')
        project = false;
    end
    if ~exist('normalize', 'var')
        normalize = 1;
    end
    if ~exist('concat', 'var')
        concat = 0;
    end
    if ~exist('novelty', 'var')
        novelty = 1;
    end
    if ~exist('which_partitions', 'var')
        which_partitions = [1,2,3];
    end
    if ~exist('which_games', 'var')
        which_games = []; % all games
    end



    [~,maskname,~] = fileparts(mask);
    %filename = sprintf('fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=%s_model=%s_%s_nsamples=100_project=%d_fast=%d_nowhiten_nofilter.mat', subj, use_smooth, glmodel, maskname, model_name, what, project, fast);
    %filename = sprintf('fit_gp_CV_HRR_cannon_repro_subj=%d_us=%d_glm=%d_mask=%s_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=%d_saveYhat=%d.mat', subj, use_smooth, glmodel, maskname, model_name, what, project, normalize, concat, novelty, fast, save_Y_hat);
    which_partitions_str = sprintf('%d', which_partitions);
    filename = sprintf('fit_gp_CV_subj=%d_us=%d_glm=%d_mask=%s_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=%d_saveYhat=%d_parts=%s.mat', subj, use_smooth, glmodel, maskname, model_name, what, project, normalize, concat, novelty, fast, save_Y_hat, which_partitions_str);
    filename = fullfile(get_mat_dir(2), filename);
    filename

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);

    %addpath(genpath('/ncf/gershman/Lab/scripts/gpml'));

    fprintf('loading BOLD for subj %d\n', subj);
    tic
    [Y, K, W, R, SPM_run_id, R_] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
    run_id = get_behavioral_run_id(subj, SPM_run_id)';
    toc

    fprintf('loading kernel for subj %d\n', subj);
    tic
    switch model_name
        case 'EMPA'
            assert(ismember(what, {'theory', 'sprite', 'interaction', 'termination', 'novelty'}));
            ker = load_HRR_kernel(subj, unique(run_id), what, normalize, concat, novelty);
        case 'DQN'
            assert(ismember(what, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2', 'all'}));
            ker = load_DQN_kernel(subj, unique(run_id), what, normalize, '');
        case 'DQN25M'
            assert(ismember(what, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2', 'all'}));
            ker = load_DQN_kernel(subj, unique(run_id), what, normalize, '25M');
        case 'PCA'
            ker = load_PCA_kernel(subj, unique(run_id), normalize);
        case 'VAE'
            ker = load_VAE_kernel(subj, unique(run_id), normalize);
        case 'game'
            ker = load_game_kernel(EXPT, subj); % GLM 1 game id features
        case 'nuisance'
            ker = load_nuisance_kernel(EXPT, subj);
        case 'state'
            ker = load_state_kernel(EXPT, subj); 
        case 'irrelevant'
            ker = load_irrelevant_kernel(EXPT, subj); 
        otherwise
            assert(false, 'invalid model name')
    end
    toc

    fprintf('Memory usage: %.3f MB\n', monitor_memory_whos);

    % whiten, filter & project out nuisance regressors
    if project
        % no white and no filter
        %Y = R_*Y;
        %ker = R_*ker*R_';
        Y = R*K*W*Y;
        ker = R*K*W*ker*W'*K'*R';
    end

    % every couple of runs form a partition
    partition_id = partition_id_from_run_id(run_id);
    assert(size(partition_id, 1) == size(Y, 1));

    % only include partitions that we care about
    partitions = unique(partition_id);
    partitions = partitions(ismember(partitions, which_partitions));
    n_partitions = length(partitions);

    disp('run_id');
    disp(run_id');
    disp('partition_id');
    disp(partition_id');
    disp('partitions');
    disp(partitions');

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

    if save_Y_hat
        Y_hat = nan(size(Y));
        Y_hat_CV = nan(size(Y));
    end

    % figure out which TRs we are including
    %
    all = ismember(partition_id, partitions); % sub select partitions
    if ~isempty(which_games)
        % sub select games
        [games, levels] = get_game_for_each_TR(subj);
        assert(length(games) == length(all));
        all = all & ismember(games, which_games);
    end
    keyboard

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
         invKi{j}] = calc_invKi(ker(all, all), x(all), s, hyp, covfun);

        ker_invKi{j} = ker(all, all) * invKi{j};

        toc

        disp('    ... for CV');
        tic

        %
        % CV
        %
        for p = 1:n_partitions % for each partition
            train = (partition_id ~= partitions(p)) & all;
            test = (partition_id == partitions(p)) & all;

            if ~any(test) || ~any(train)
                % empty partition 
                disp(['kernel skipping empty partition', num2str(partitions(p))]);
                continue;
            end

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
        all = ismember(partition_id, partitions); % all trials
        train = all;
        test = all; % test data = train data here

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
            train = (partition_id ~= partitions(p)) & all;
            test = (partition_id == partitions(p)) & all;

            if ~any(test) || ~any(train)
                % empty partition 
                disp(['skipping empty partition', num2str(partitions(p))]);
                continue;
            end

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

        if save_Y_hat
            Y_hat(:,i) = y_hat;
            Y_hat_CV(:,i) = y_hat_CV;
        end

        %figure;
        %hold on;
        %plot(y);
        %plot(y_hat);
        %plot(y_hat_CV);
        %legend({'y', 'y_hat', 'y_hat_CV'}, 'interpreter', 'none');

        %figure;
        %hold on;
        %plot(y);
        %plot(ceil_y_hat_CV);
        %legend({'y', 'ceil_y_hat', 'ceil_y_hat_CV'}, 'interpreter', 'none');

        %break
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
    if save_Y_hat
        save(filename, 'sigma', 'logmarglik', 'logpredlik', 'R2', 'adjR2', 'r', 'MSE', 'SMSE', ... 
                       'sigma_CV', 'logmarglik_CV', 'logpredlik_CV', 'R2_CV', 'adjR2_CV', 'r_CV', 'MSE_CV', 'SMSE_CV', ...
                       'ceil_sigma', 'ceil_logmarglik', 'ceil_logpseudolik', 'ceil_R2', 'ceil_adjR2', 'ceil_r', 'ceil_MSE', 'ceil_SMSE', ...
                       'ceil_sigma_CV', 'ceil_logmarglik_CV', 'ceil_logpseudolik_CV', 'ceil_R2_CV', 'ceil_adjR2_CV', 'ceil_r_CV', 'ceil_MSE_CV', 'ceil_SMSE_CV', ...
                       'subj', 'use_smooth', 'glmodel', 'mask', 'what', 'sigmas', 'n', 'partition_id', 'partitions', 'Y', 'Y_hat_CV', 'Y_hat', ...
        '-v7.3');
    else
        save(filename, 'sigma', 'logmarglik', 'logpredlik', 'R2', 'adjR2', 'r', 'MSE', 'SMSE', ... 
                       'sigma_CV', 'logmarglik_CV', 'logpredlik_CV', 'R2_CV', 'adjR2_CV', 'r_CV', 'MSE_CV', 'SMSE_CV', ...
                       'ceil_sigma', 'ceil_logmarglik', 'ceil_logpseudolik', 'ceil_R2', 'ceil_adjR2', 'ceil_r', 'ceil_MSE', 'ceil_SMSE', ...
                       'ceil_sigma_CV', 'ceil_logmarglik_CV', 'ceil_logpseudolik_CV', 'ceil_R2_CV', 'ceil_adjR2_CV', 'ceil_r_CV', 'ceil_MSE_CV', 'ceil_SMSE_CV', ...
                       'subj', 'use_smooth', 'glmodel', 'mask', 'what', 'sigmas', 'n', 'partition_id', 'partitions', ...
        '-v7.3');
    end

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


