% aggregate results from fit_gp_CV.m
% copy of agg_gp.m

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

    %filename = sprintf('mat/fit_gp_HRR_subj=%d_us=%d_glm=21_mask=mask_%s.mat', subj_id, use_smooth, what);
    filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=21_mask=mask_%s.mat', subj_id, use_smooth, what);
    load(filename, 'n', 'logmarglik', 'null_logmarglik', 'adjR2', 'R2_CV', 'null_adjR2', 'null_R2_CV', 'sigma', 'mask', 'r', 'null_r', 'r_CV', 'null_r_CV', 'mask');

    % calc stuff
    k = 1; % 1 param = sigma
    n = 1698; % # TRs
    bic = k*log(n) - 2*logmarglik;
    null_bic = k*log(n) - 2*null_logmarglik;

    if s == 1
        logBF = nan(length(subjects), size(r,2));
        rs = nan(length(subjects), size(r,2));
        log_group_marglik = nan(1, size(r,2));
        null_log_group_marglik = nan(1, size(r,2));
    end

    log_group_marglik = log_group_marglik + logmarglik;
    null_log_group_marglik = null_log_group_marglik + null_logmarglik;

    % log Bayes factor
    logBF(s,:) = logmarglik - null_logmarglik;

    % CV Pearson r's across subjects
    rs(s,:) = mean(r_CV, 1);
    null_rs(s,:) = mean(null_r_CV, 1);

    % Fisher z transform
    zs = atanh(rs);
    null_zs = atand(null_rs);
end

% group log BF = sum of log BF across subjects (Stephan et al. 2009)
logGBF = sum(logBF,1);

% t-test across subjects
[h,p,ci,stats] = ttest(zs, null_zs);
ts = stats.tstat;

% create SPMs
%
logGBFmap = zeros(size(mask));
logGBFmap(mask) = logGBF;

tmap = zeros(size(mask));
tmap(mask) = ts;

loglikmap = zeros(size(mask));
loglikmap(mask) = log_group_marglik;

null_loglikmap = zeros(size(mask));
null_loglikmap(mask) = null_log_group_marglik;


filename = sprintf('mat/agg_gp_CV_us=%d_%s.mat', use_smooth, what);
filename
save(filename);

