% aggregate results from fit_gp_CV.m
% copy of agg_gp.m

use_smooth = true;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);
%subjects = 1:2:32; % odd

model_name = 'nuisance';
what = '';
project = 0;
glmodel = 1;

%agg_filename = sprintf('mat/agg_gp_CV_us=%d_glm=%d_%s_fast.mat', use_smooth, glmodel, what);
%agg_filename = sprintf('mat/agg_gp_CV_us=%d_glm=%d_model=EMPA_%s_nsamples=100_fast_WTF.mat', use_smooth, glmodel, what);
%agg_filename = sprintf('mat/agg_gp_CV_us=%d_glm=%d_model=EMPA_%s_nsamples=100_fast.mat', use_smooth, glmodel, what);
agg_filename = sprintf('mat/agg_gp_CV_us=%d_glm=%d_model=%s_%s_nsamples=100_project=%d_fast=1.mat', use_smooth, glmodel, model_name, what, project);
agg_filename

for s = 1:length(subjects)
    subj_id = subjects(s);

    %filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_%s_fast.mat', subj_id, use_smooth, glmodel, what);
    %filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=EMPA_%s_nsamples=100_fast.mat', subj_id, use_smooth, glmodel, what);
    %filename = sprintf('/Volumes/fMRI-2/Mac_mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=EMPA_%s_nsamples=100_fast.mat', subj_id, use_smooth, glmodel, what);
    %filename = sprintf('/Volumes/fMRI-2/Mac_mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=EMPA_%s_nsamples=100_project=%d_fast=1_.mat', subj_id, use_smooth, glmodel, what, project);
    filename = sprintf('/Volumes/fMRI-2/Mac_mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_nsamples=100_project=%d_fast=1.mat', subj_id, use_smooth, glmodel, model_name, what, project);
    filename
    %filename = sprintf('mat/fit_gp_CV_noRKW_HRR_subj=%d_us=%d_glm=21_mask=mask_%s.mat', subj_id, use_smooth, what);

    load(filename, 'n', 'logmarglik', 'adjR2', 'R2_CV', 'sigma', 'mask', 'r', 'r_CV', 'mask', 'MSE', 'SMSE', 'MSE_CV', 'SMSE_CV', 'partition_id');

    % calc stuff
    k = 1; % 1 param = sigma
    n = length(partition_id); % # TRs
    bic = k*log(n) - 2*logmarglik;

    if s == 1
        logBF = nan(length(subjects), size(r,2));

        rs = nan(length(subjects), size(r,2));

        adjR2s = nan(length(subjects), size(r,2));

        R2s = nan(length(subjects), size(r,2));

        MSEs = nan(length(subjects), size(r,2));
        SMSEs = nan(length(subjects), size(r,2));

        log_group_marglik = zeros(1, size(r,2));
    end

    log_group_marglik = log_group_marglik + logmarglik;

    % log Bayes factor
    %logBF(s,:) = logmarglik - null_logmarglik;

    % MSEs
    MSEs(s,:) = mean(MSE_CV, 1);
    SMSEs(s,:) = mean(SMSE_CV, 1);

    % R2s
    R2s(s,:) = mean(R2_CV, 1);

    adjR2s(s,:) = adjR2;
e
    % CV Pearson r's across subjects
    rs(s,:) = mean(r_CV, 1);

end

% Fisher z transform
zs = atanh(rs);

% group log BF = sum of log BF across subjects (Stephan et al. 2009)
%logGBF = sum(logBF,1);

% t-test Pearson corr across subjects
[h,p,ci,stats] = ttest(zs);
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
%logGBFmap = zeros(size(mask));
%logGBFmap(mask) = logGBF;

R2map = zeros(size(mask));
R2map(mask) = mean(R2s,1);

adjR2map = zeros(size(mask));
adjR2map(mask) = mean(adjR2s,1);


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


%{
Wmap = zeros(size(mask));
Ws(ps < 0.01) = nan;
Wmap(mask) = Ws;
%}


agg_filename
save(agg_filename);

