function [se, m] = wse(X)
    
    % Within-subject error, following method of Cousineau (2005).
    %
    % USAGE: [se, m] = wse(X)
    %
    % INPUTS:
    %   X - [N x D] data with N subjects and D observations
    %
    % OUTPUTS:
    %   se - [1 x D] within-subject standard errors
    %   m - [1 x D] means
    %
    % Sam Gershman, June 2015
    
    m = squeeze(nanmean(X));
    X = bsxfun(@minus,X,nanmean(X,2));
    N = sum(~isnan(X));
    se = bsxfun(@rdivide,nanstd(X),sqrt(N));
    se = squeeze(se);
