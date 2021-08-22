function stars = significance(p)
    % p-value to *'s

    if p > 0.05
        stars = 'ns';
    elseif p > 0.01
        stars = '*';
    else
        stars = ['*', repmat('*', 1, floor(-log10(p)) - 1)];
    end
end

