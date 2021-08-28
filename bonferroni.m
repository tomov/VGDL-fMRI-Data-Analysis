function ps_corr = bonferroni(ps)
    % apply Bonferroni correction to a set of P values
     ps_corr = 1 - (1 - ps) .^ length(ps);
