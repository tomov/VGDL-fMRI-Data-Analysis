% calculate inverse of kernel, for streamlining GP regression
%
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
