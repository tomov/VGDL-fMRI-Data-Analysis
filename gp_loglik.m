function [loglik, grad] = gp_loglik(K, y, sigma, neg)

    % GP regression marginal likelihood
    % [loglik, grad] = gp_loglik(K, y, sigma)
    %
    % K - [n x n] K(X,X) kernel = X * Sigma_p * X'
    % y - [n x 1] training output
    % sigma - noise std dev
    % neg - whether to flip the signs (for fmincon)

    if ~exist('neg', 'var')
        neg = false;
    end

    I = eye(size(K));
    Ki = K + sigma^2 * I;
    invKi = Ki^(-1);

    loglik = -0.5 * y' * invKi * y - 0.5 * log(det(Ki)) - length(y)/2 * log(2 * pi); % Eq. 2.30 from Rasmussen
    if neg
        loglik = -loglik;
    end

    if nargout > 1
        grad = 0.5 * y' * invKi * 2*sigma*I * invKi * y - 0.5 * trace(invKi * 2*sigma*I); % Eq 5.9
        if neg
            grad = -grad;
        end
    end
