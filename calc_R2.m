function [R2, adjR2] = calc_R2(y, y_hat, p)
    % calculate R^2 https://en.wikipedia.org/wiki/Coefficient_of_determination
    %
    % p = # of free parameters (e.g. sigma)
    %
    n = length(y);

    SStot = sum((y - mean(y)).^2); 
    assert(immse(SStot, var(y) * (n - 1)) < 1e-15);
    SSres = sum((y - y_hat).^2);

    R2 = 1 - SSres/SStot;

    adjR2 = 1 - (1 - R2) * (n - 1) / (n - p - 1);
end


