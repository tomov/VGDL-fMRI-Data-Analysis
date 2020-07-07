% aggregate results from fit_gp.m

use_smooth = true;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);
what = 'theory';


for s = 1:length(subjects)
    subj_id = subjects(s);

    filename = sprintf('mat/fit_gp_HRR_subj=%d_us=%d_glm=21_mask=mask_%s.mat', subj_id, use_smooth, what);
    load(filename, 'n', 'loglik', 'adjR2', 'sigma', 'mask');

    % calc stuff
    k = 1; % 1 param = sigma
    bic = k*log(n) - 2*loglik;

    % aggregate
    if s == 1
        mean_adjR2 = adjR2;
        group_bic = bic;
        mean_sigma = sigma;
    else
        mean_adjR2 = mean_adjR2 + adjR2;
        group_bic = group_bic + bic;
        mean_sigma = mean_sigma + sigma;
    end
end

mean_adjR2 = mean_adjR2 / length(subjects);
mean_sigma = mean_sigma / length(subjects);

filename = sprintf('mat/agg_gp_us=%d_%s.mat', use_smooth, what);
filename
save(filename);

