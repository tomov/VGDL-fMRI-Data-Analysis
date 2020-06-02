rng(334);

X = [ones(50,1) rand(50,3)];
b = [0; rand(3,1)];
e = normrnd(0, 0.01, size(X,1), 1);
y = X * b + e;

% test OLS

b0 = glmfit(X,y, 'normal','constant', 'off');

beta = (X' * X)^(-1) * X' * y;  % OLS
assert(immse(beta, b0) < 1e-5);

% test ridge with k = 0 i.e. OLS

b2 = ridge(y,X(:,2:end),0,0);
assert(immse(b2, b0) < 1e-5);

% test ridge

k = 0.001;
b2 = ridge(y,X(:,2:end),k,0);

%beta2 = (X' * X + k * eye(size(X,2)))^(-1) * X' * y; % doesn't work; see https://www.mathworks.com/matlabcentral/answers/29373-ridge-regression-coefficient-question
beta2 = (X' * X + k * var(diag(X)))^(-1) * X' * y; % note that ridge regression is not scale invariant, that is, X(:,i) * 2 doesn't make b(i) * 2
assert(immse(b2, beta) < 1e-5);
