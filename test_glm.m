clear all;

addpath(genpath('/Users/momchil/Dropbox/Research/libs/NearestSymmetricPositiveDefinite/')); 

rng(334);

sigma = 0.01;
tau = 3;
X = [ones(50,1) rand(50,3)];
b = [0; normrnd(0, tau, 3, 1)];
e = normrnd(0, sigma, size(X,1), 1);
y = X * b + e;

Xnew = [ones(50,1) rand(50,3)];
ynew = Xnew * b + normrnd(0, sigma, size(X,1), 1);

%% test OLS

% fit glm and manually
beta = (X' * X)^(-1) * X' * y;  % OLS
b0 = glmfit(X,y, 'normal','constant', 'off');

assert(immse(beta, b0) < 1e-10);

% predict
Ynew_0 = X * beta;
Ynew0 = X * b0;

assert(immse(Ynew0, Ynew_0) < 1e-10);


%% test ridge with k = 0 i.e. OLS

% fit ridge
b1 = ridge(y,X(:,2:end),0,0); % automatically includes intercept

assert(immse(b1, b0) < 1e-10);

% predict
Ynew1 = X * b1;

assert(immse(Ynew1, Ynew0) < 1e-5);


%% test ridge

k = 0.001;

% fit ridge and manually
beta2 = (X' * X + k * var(diag(X)))^(-1) * X' * y; % note that ridge regression is not scale invariant, that is, X(:,i) * 2 doesn't make b(i) * 2
%beta2 = (X' * X + k * eye(size(X,2)))^(-1) * X' * y; % doesn't work; see https://www.mathworks.com/matlabcentral/answers/29373-ridge-regression-coefficient-question
b2 = ridge(y,X(:,2:end),k,0);

assert(immse(b2, beta) < 1e-5);

% predict
Ynew_2 = Xnew * beta2;
Ynew2 = Xnew * b2;

assert(immse(Ynew2, Ynew_2) < 1e-5);

%% test GP

% fit classic ridge first
lambda = sigma^2 / tau^2; % in accordance with Bayesian interpretation; see https://statisticaloddsandends.wordpress.com/2018/12/29/bayesian-interpretation-of-ridge-regression/
beta3 = (X' * X + lambda * eye(size(X,2)))^(-1) * X' * y; 

% fit ridge according to Rasmussen (2006) eq. 2.8
% notice X is transposed compared to their notation
Sigma_p = tau.^2 * eye(size(X,2));
A = 1/sigma^2 * X' * X + Sigma_p^(-1);
beta4 = 1/sigma^2 * A^(-1) * X' * y;

assert(immse(beta3, beta4) < 1e-15); % exactly identical

% predict
Ynew3 = Xnew * beta3;
Ynew4 = Xnew * beta4;

assert(immse(Ynew3, Ynew4) < 1e-10);

% predict using kernels, Rasmussen (2006) eq. 2.12; also see Eq 2.25 and 2.26
K = X * Sigma_p * X';
I = eye(size(X,1));
K = nearestSPD(K);  % find nearest symmetric positive definite matrix (it's not b/c of numerical issues, floating points, etc.)
Ynew5 = Xnew * Sigma_p * X' * (K + sigma^2 * I) ^ (-1) * y;
varYnew5 = Xnew * Sigma_p * Xnew' - Xnew * Sigma_p * X' * (K + sigma^2 * I) ^ (-1) * X * Sigma_p * Xnew';

assert(immse(Ynew5, Ynew3) < 1e-5);

% calculate log marginal likelihood (eq 2.30) explicitly and as the MVN pdf
margLik = (2 * pi)^(-length(y)/2) * det(K + sigma^2 * I)^(-1/2) * exp(-0.5 * y' * (K + sigma^2 * I)^(-1) * y);
logMargLik = -0.5 * y' * (K + sigma^2 * I)^(-1) * y - 0.5 * log(det(K + sigma^2 * I)) - length(y)/2 * log(2 * pi);

margLik2 = mvnpdf(y', zeros(size(y')), K + sigma^2 * I);
logMargLik2 = log(margLik2);

assert(abs(1 - margLik / margLik2) < 1e-5); % look at ratio; b/c they're huge
assert(immse(logMargLik, logMargLik2) < 1e-10); 

% fit GP 
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
[Ynew5 ~] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y, xs);

assert(immse(Ynew5, Ynew3) < 1e-10); % for some reason, doesn't work..........

close all; plot(Ynew3); hold on; plot(Ynew5);
