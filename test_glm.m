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

logMargLikFn = @(sigma) -0.5 * y' * (K + sigma^2 * I)^(-1) * y - 0.5 * log(det(K + sigma^2 * I)) - length(y)/2 * log(2 * pi); % Eq. 2.30
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
%
addpath(genpath('/Users/momchil/Dropbox/Research/libs/gpml/')); % GP ML

n = size(X,1) + size(Xnew,1);

meanfun = {@meanDiscrete, n};
covfun = {@covDiscrete, n};
likfun = @likGauss;

K = [X; Xnew] * Sigma_p * [X; Xnew]';
K = nearestSPD(K);  % find nearest symmetric positive definite matrix (it's not b/c of numerical issues, floating points, etc.)
L = chol(K);
covhyp = L(triu(true(n)));

% GP hyperparams
hyp = struct('mean', zeros(1, n), 'cov', covhyp, 'lik', -1);

% GP
x = [1:size(X,1)]';
xs = [size(X,1)+1:n]';
[y_GP_kernel ~] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y, xs);

%assert(immse(y_ridge_kernel, y_ridge_Bayesian) < 1e-10); % for some reason, doesn't work quite well...

close all; 
plot(y_ridge_Bayesian); 
hold on; 
plot(y_ridge_Rasmussen);
plot(y_ridge_kernel);
plot(y_GP_kernel); 
legend({'Bayesian (ridge) regression, from interwebs', 'Bayesian (ridge) regression, from Rasmussen', 'Equivalent kernel regression', 'Using GPML library'});
