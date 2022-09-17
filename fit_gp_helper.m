% helper for fit_gp_CV.m
%
% fit gp sigmas manually; optionally sanity check with gp() library
%
function [sigma, logmarglik, logpredlik, y_hat, R2, adjR2, r, MSE, SMSE] = fit_gp_helper(x, y, train, test, ker, hyp, meanfun, covfun, likfun, sigmas, invKi, ldB2, sn2, solveKiW, ker_invKi, fast, debug)

    n = sum(train);

    % for each sigma, compute marginal likelihood = loglik = -NLZ
    %
    for j = 1:length(sigmas)
        s = sigmas(j);
            
        % from infGaussLik.m
        keyboard
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
    if fast
        logpredlik = nan; % it's slow
    else
        hyp.lik = log(sigma);
        [~,~,~,~,lp] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x(train), y(train), x(test), y(test));
        logpredlik = sum(lp);
    end

    % posterior predictive on test dataset (Eq. 2.23 in Rasmussen)
    y_hat = ker_invKi{j} * y(train);

    % R^2
    %assert(all(train == test) || ~any(train == test)); % <-- not true if using some partitions
    if all(train == test)
        p = 1; % # 1 (hyper)param = sigma; b/c we fit on same data 
    else
        p = 0; % 0 params, b/c evaluating on held out data
    end
    [R2, adjR2] = calc_R2(y(test), y_hat, 1); % TODO p here?

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


