function [l, g] = test_fit_GP(sigma, K, y)
    % compute LME and its gradient for GP
    % see test_glm.m

    I = eye(size(K));
    Ky = K + sigma^2 * I;
    invKy = Ky^(-1);

    if nargout == 1
        l = gp_loglik(K, y, sigma);
        l = -l;
    else
        assert(nargout == 2);
        [l, g] = gp_loglik(K, y, sigma);
        l = -l;
        g = -g;
    end
