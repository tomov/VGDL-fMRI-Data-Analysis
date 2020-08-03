% aggregate results from fit_gp_CV.m
% copy of agg_gp.m

use_smooth = true;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);
what = 'termination';
glmodel = 9;

agg_filename = sprintf('mat/agg_gp_CV_us=%d_glm=%d_%s.mat', use_smooth, glmodel, what);

for s = 1:length(subjects)
    subj_id = subjects(s);

    filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_%s.mat', subj_id, use_smooth, glmodel, what);
    %filename = sprintf('mat/fit_gp_CV_noRKW_HRR_subj=%d_us=%d_glm=21_mask=mask_%s.mat', subj_id, use_smooth, what);
    load(filename, 'n', 'logmarglik', 'null_logmarglik', 'adjR2', 'R2_CV', 'null_adjR2', 'null_R2_CV', 'sigma', 'mask', 'r', 'null_r', 'r_CV', 'null_r_CV', 'mask', 'MSE', 'null_MSE', 'SMSE', 'null_SMSE', 'MSE_CV', 'SMSE_CV', 'null_MSE_CV', 'null_SMSE_CV');

    % calc stuff
    k = 1; % 1 param = sigma
    n = 1698; % # TRs
    bic = k*log(n) - 2*logmarglik;
    null_bic = k*log(n) - 2*null_logmarglik;

    if s == 1
        logBF = nan(length(subjects), size(r,2));

        rs = nan(length(subjects), size(r,2));
        null_rs = nan(length(subjects), size(r,2));

        adjR2s = nan(length(subjects), size(r,2));

        R2s = nan(length(subjects), size(r,2));
        null_R2s = nan(length(subjects), size(r,2));

        MSEs = nan(length(subjects), size(r,2));
        SMSEs = nan(length(subjects), size(r,2));
        null_MSEs = nan(length(subjects), size(r,2));
        null_SMSEs = nan(length(subjects), size(r,2));

        log_group_marglik = nan(1, size(r,2));
        null_log_group_marglik = nan(1, size(r,2));
    end

    log_group_marglik = log_group_marglik + logmarglik;
    null_log_group_marglik = null_log_group_marglik + null_logmarglik;

    % log Bayes factor
    logBF(s,:) = logmarglik - null_logmarglik;

    % MSEs
    MSEs(s,:) = mean(MSE_CV, 1);
    SMSEs(s,:) = mean(SMSE_CV, 1);
    null_MSEs(s,:) = mean(null_MSE_CV, 1);
    null_SMSEs(s,:) = mean(null_SMSE_CV, 1);

    % R2s
    R2s(s,:) = mean(R2_CV, 1);
    null_R2s(s,:) = mean(null_R2_CV, 1);

    adjR2s(s,:) = adjR2;

    % CV Pearson r's across subjects
    rs(s,:) = mean(r_CV, 1);
    null_rs(s,:) = mean(null_r_CV, 1);

end

% Fisher z transform
zs = atanh(rs);
null_zs = atand(null_rs);

% group log BF = sum of log BF across subjects (Stephan et al. 2009)
logGBF = sum(logBF,1);

% t-test Pearson corr across subjects
[h,p,ci,stats] = ttest(zs);
%[h,p,ci,stats] = ttest(zs, null_zs);
ts = stats.tstat;

%{
% don't: not Gaussian

% t-test SMSEs across subjects
[h,p,ci,stats] = ttest(SMSEs, null_SMSEs);
SMSE_ts = stats.tstat;

% t-test MSEs across subjects
[h,p,ci,stats] = ttest(MSEs, null_MSEs);
MSE_ts = stats.tstat;

% t-test R2 across subjects
[h,p,ci,stats] = ttest(R2s, null_R2s);
R2_ts = stats.tstat;
%}

% sign-rank test MSEs across subjects 
% not very informative...
%{
tic
ps = nan(size(ts));
Ws = nan(size(ts));
for j = 1:size(MSEs,2)
    [h,ps(j),stat] = signrank(MSEs(:,j), null_MSEs(:,j));
    Ws(j) = stat.signedrank;
end
toc
%}

% create SPMs
%
logGBFmap = zeros(size(mask));
logGBFmap(mask) = logGBF;

R2map = zeros(size(mask));
R2map(mask) = mean(R2s,1);

adjR2map = zeros(size(mask));
adjR2map(mask) = mean(adjR2s,1);

null_R2map = zeros(size(mask));
null_R2map(mask) = mean(null_R2s,1);

tmap = zeros(size(mask));
tmap(mask) = ts;

%{
SMSE_tmap = zeros(size(mask));
SMSE_tmap(mask) = SMSE_ts;

MSE_tmap = zeros(size(mask));
MSE_tmap(mask) = MSE_ts;

R2_tmap = zeros(size(mask));
R2_tmap(mask) = R2_ts;
%}

loglikmap = zeros(size(mask));
loglikmap(mask) = log_group_marglik;

null_loglikmap = zeros(size(mask));
null_loglikmap(mask) = null_log_group_marglik;

%{
Wmap = zeros(size(mask));
Ws(ps < 0.01) = nan;
Wmap(mask) = Ws;
%}


agg_filename
save(agg_filename);

