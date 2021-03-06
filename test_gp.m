clear all;

% Gaussian process vs. ridge regression vs. OLS
% GP

addpath(genpath('/Users/momchil/Dropbox/Research/libs/NearestSymmetricPositiveDefinite/')); 

rng(334);

sigma = 0.01;
tau = 3;
X = [ones(50,1) rand(50,3)];
b = [0; normrnd(0, tau, 3, 1)];
e = normrnd(0, sigma, size(X,1), 1);
y = X * b + e;

Xnew = [ones(30,1) rand(30,3)];
ynew = Xnew * b + normrnd(0, sigma, size(Xnew,1), 1);

% ---------------------------------------------
%               test OLS
% ---------------------------------------------

% fit glm and manually
beta_OLS_manual = (X' * X)^(-1) * X' * y;  % OLS
beta_OLS_glmfit = glmfit(X,y, 'normal','constant', 'off');

assert(immse(beta_OLS_manual, beta_OLS_glmfit) < 1e-10);

% predict
y_OLS_manual = X * beta_OLS_manual;
y_OLS_glmfit = X * beta_OLS_glmfit;

assert(immse(y_OLS_glmfit, y_OLS_manual) < 1e-10);


% ---------------------------------------------
%        test ridge with k = 0 i.e. OLS
% ---------------------------------------------

% fit ridge
beta_OLS_ridge = ridge(y,X(:,2:end),0,0); % automatically includes intercept

assert(immse(beta_OLS_ridge, beta_OLS_manual) < 1e-10);

% predict
y_OLS_ridge = X * beta_OLS_ridge;

assert(immse(y_OLS_ridge, y_OLS_glmfit) < 1e-5);


% ---------------------------------------------
%                   test ridge
% ---------------------------------------------

k = 0.001;

% fit ridge and manually
beta_ridge_manual = (X' * X + k * var(diag(X)))^(-1) * X' * y; % note that ridge regression is not scale invariant, that is, X(:,i) * 2 doesn't make b(i) * 2
%beta_ridge_manual = (X' * X + k * eye(size(X,2)))^(-1) * X' * y; % doesn't work; see https://www.mathworks.com/matlabcentral/answers/29373-ridge-regression-coefficient-question
beta_ridge = ridge(y,X(:,2:end),k,0);

assert(immse(beta_ridge, beta_ridge_manual) < 1e-5);

% predict
y_ridge_manual = Xnew * beta_ridge_manual;
y_ridge = Xnew * beta_ridge;

assert(immse(y_ridge, y_ridge_manual) < 1e-5);

% ---------------------------------------------
%                     test GP
% ---------------------------------------------

%% fit classic ridge first
%
lambda = sigma^2 / tau^2; % in accordance with Bayesian interpretation; see https://statisticaloddsandends.wordpress.com/2018/12/29/bayesian-interpretation-of-ridge-regression/
beta_ridge_Bayesian = (X' * X + lambda * eye(size(X,2)))^(-1) * X' * y; 

%% fit ridge according to Rasmussen (2006) eq. 2.8
%
% notice X is transposed compared to their notation
Sigma_p = tau.^2 * eye(size(X,2));
A = 1/sigma^2 * X' * X + Sigma_p^(-1);
beta_ridge_Rasmussen = 1/sigma^2 * A^(-1) * X' * y; 

assert(immse(beta_ridge_Bayesian, beta_ridge_Rasmussen) < 1e-15); % exactly identical

%% predict ridge
%
y_ridge_Bayesian = Xnew * beta_ridge_Bayesian;
y_ridge_Rasmussen = Xnew * beta_ridge_Rasmussen;

assert(immse(y_ridge_Bayesian, y_ridge_Rasmussen) < 1e-10);



%% predict using kernels, Rasmussen (2006) eq. 2.12; also see Eq 2.25 and 2.26
%
K = X * Sigma_p * X';
K_new = X * Sigma_p * Xnew';
K_new_new = Xnew * Sigma_p * Xnew';
K = nearestSPD(K);  % find nearest symmetric positive definite matrix (it's not b/c of numerical issues, floating points, etc.)
I = eye(size(X,1));

predFn = @(sigma) Xnew * Sigma_p * X' * (K + sigma^2 * I) ^ (-1) * y; % Eq. 2.12, Eq. 2.25

y_ridge_kernel_1 = predFn(sigma); % mean prediction
vary_ridge_kernel_1 = Xnew * Sigma_p * Xnew' - Xnew * Sigma_p * X' * (K + sigma^2 * I) ^ (-1) * X * Sigma_p * Xnew'; % variance of prediction; Eq. 2.12, Eq. 2.26


[y_ridge_kernel, vary_ridge_kernel] = gp_pred(K, y, K_new, K_new_new, sigma);
assert(immse(y_ridge_kernel, y_ridge_kernel_1) < 1e-15);
assert(immse(vary_ridge_kernel, vary_ridge_kernel_1) < 1e-15);

assert(immse(y_ridge_kernel, y_ridge_Bayesian) < 1e-5);

logMargLikFn = @(sigma) -0.5 * y' * (K + sigma^2 * I)^(-1) * y - 0.5 * log(det(K + sigma^2 * I)) - length(y)/2 * log(2 * pi); % Eq. 2.30 TODO is there a sigma^2 after the 2 * pi? not in the book, but yes in infGaussLik
assert(immse(logMargLikFn(sigma), -test_fit_GP(sigma, K, y)) < 1e-10);



% calculate log marginal likelihood (eq 2.30) manually and as the MVN pdf
%
margLik = (2 * pi)^(-length(y)/2) * det(K + sigma^2 * I)^(-1/2) * exp(-0.5 * y' * (K + sigma^2 * I)^(-1) * y);
logMargLik = logMargLikFn(sigma)

margLik2 = mvnpdf(y', zeros(size(y')), K + sigma^2 * I);
logMargLik2 = log(margLik2);

logMargLik3 = gp_loglik(K, y, sigma);

assert(abs(1 - margLik / margLik2) < 1e-5); % look at ratio; b/c they're huge
assert(immse(logMargLik, logMargLik2) < 1e-10); 
assert(immse(logMargLik3, logMargLik2) < 1e-10); 

% fit GP hyperparams manually 

options = optimoptions('fmincon','SpecifyObjectiveGradient',true)
sigma_hat = fmincon(@(sigma) gp_loglik(K, y, sigma, true), 1, -1, 0, [], [], [], [], [], options);
sigma_hat
sigma

y_ridge_kernel_hyperparam = predFn(sigma_hat);

assert(immse(y_ridge_kernel_hyperparam, y_ridge_Bayesian) < 1e-5);

%% fit GP using Rasmussen's library
% http://www.gaussianprocess.org/gpml/code/matlab/doc/
%

% init stuff
%
n = size(X,1) + size(Xnew,1);

meanfun = {@meanDiscrete, n};
covfun = {@covDiscrete, n};
likfun = @likGauss;

K_all = [X; Xnew] * Sigma_p * [X; Xnew]';
K_all = nearestSPD(K_all);  % find nearest symmetric positive definite matrix (it's not b/c of numerical issues, floating points, etc.)
L = cholcov(K_all); % not chol b/c sometimes positive semi-definite
L(1:(n+1):end) = log(diag(L));  % see covDiscrete.m
covhyp = L(triu(true(n)));  % see covDiscrete.m

% GP hyperparams
hyp = struct('mean', zeros(1, n), 'cov', covhyp, 'lik', log(sigma));


% GP prediction
x = [1:size(X,1)]';
xnew = [size(X,1)+1:n]';
[y_GP_kernel, ~, ~, ~, lp] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y, xnew, ynew);
logPred_GP_kernel = sum(lp);
logPred_GP_kernel

% sanity
[K_gp,dK] = feval(@covDiscrete,n,hyp.cov,x);
assert(immse(K_gp, K) < 1e-15); % should be identical

% from infLikGauss.m
Kf_gp = apx(hyp, covfun, x, []);
assert(immse(Kf_gp.mvm(eye(size(K))), K) < 1e-15); % should be identical

assert(immse(y_ridge_kernel, y_ridge_Bayesian) < 1e-15);

% GP marginal log likelihood
negLogMargLik_GP = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);
logMargLik_GP = -negLogMargLik_GP;
assert(immse(logMargLik_GP, logMargLik) < 1e-15);

%% minimize GP
% http://www.gaussianprocess.org/gpml/code/matlab/doc/

% fit all hyperparams -- not good; we don't want to fit kernel & mean fun
%hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfun, covfun, likfun, x, y);

% ...instead, fit noise only (TODO too slow!)
for i = 1:n
    prior.mean{i} = {@priorDelta};
end
for i = 1:length(covhyp)
    prior.cov{i} = {@priorDelta};
end
inf = {@infPrior, @infGaussLik, prior};
hyp.lik = log(10); % confuse it by starting off the true log(sigma)
hyp2 = minimize(hyp, @gp, -100, inf, meanfun, covfun, likfun, x, y);
% these should be unchanged 
assert(immse(hyp2.cov, hyp.cov) < 1e-20);
assert(immse(hyp2.mean, hyp.mean) < 1e-20);

fprintf('sigma = %.4f vs. inferred sigma = %.4f\n', sigma, exp(hyp2.lik));

hyp2
nlz = gp(hyp2, @infGaussLik, meanfun, covfun, likfun, x, y);
nlz

[y_GP_kernel_fit, ~, ~, ~, lp] = gp(hyp2, @infGaussLik, meanfun, covfun, likfun, x, y, xnew, ynew);
logPred_GP_kernel_fit = sum(lp);
logPred_GP_kernel_fit


%% plot stuff
close all; 
hold on; 
plot(ynew); 
plot(y_ridge_Bayesian + rand(size(ynew)) * 0.1); 
plot(y_ridge_Rasmussen + rand(size(ynew)) * 0.1);
plot(y_ridge_kernel + rand(size(ynew)) * 0.1);
plot(y_GP_kernel + rand(size(ynew)) * 0.1); 
plot(y_GP_kernel_fit + rand(size(ynew)) * 0.1); 
legend({'truth', 'Bayesian (ridge) regression, from interwebs', 'Bayesian (ridge) regression, from Rasmussen', 'Equivalent kernel regression', 'Using GPML library', 'GP fit hyperparams'});

%}


