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


