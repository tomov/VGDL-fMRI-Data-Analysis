function [l, g] = test_fit_GP(sigma, K, y)
    % see test_glm.m

    I = eye(size(K));
    Ky = K + sigma^2 * I;
    invKy = Ky^(-1);

    logMargLikFn = @(sigma) -0.5 * y' * invKy * y - 0.5 * log(det(Ky)) - length(y)/2 * log(2 * pi); % Eq. 2.30 from Rasmussen
    l = -logMargLikFn(sigma);

    if nargout > 1
        logMargLikFn_grad = @(sigma) 0.5 * y' * invKy * 2*sigma*I * invKy * y - 0.5 * trace(invKy * 2*sigma*I); % Eq 5.9
        g = -logMargLikFn_grad(sigma);
    end
