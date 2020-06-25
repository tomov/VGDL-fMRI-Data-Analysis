function [y_mean, y_var] = gp_pred(K, y, K_new, K_new_new, sigma)

    % GP regression prediction
    % y_new = gp_pred(K, y, K_new, K_new_new, sigma)
    %
    % K - [n x n] K(X,X) kernel = X * Sigma_p * X'
    % y - [n x 1] training output
    % K_new - [n x k] K(X,X_new) kernel = X * Sigma_p * X_new'
    % K_new_new - [k x k] K(X_new,X_new) kernel = X_new * Sigma_p * X_new'
    % sigma - noise std dev

    % Eq. 2.23; also see Eq. 2.12, Eq. 2.25 from Rasmussen book
    %
    I = eye(size(K));
    K = K + sigma^2 * I;
    invK = Ki^(-1);

    invKi

    y_mean = K_new' * invKi * y; % Eq. 2.23; also Eq. 2.12, Eq. 2.25
    y_var = K_new_new - K_new' * invKi * K_new; % variance of prediction; Eq. 2.24; also Eq. 2.12, Eq. 2.26
