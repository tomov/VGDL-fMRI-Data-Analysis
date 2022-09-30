% aggregate results from fit_gp_CV.m
% copy of agg_gp.m
close all;
clear all;


use_smooth = true;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);
%subjects = 1:2:32; % odd

%model_name = 'PCA';
%model_name = 'state';
%model_name = 'irrelevant';
%model_name = 'DQN';
%model_name = 'DQN25M';
model_name = 'DQN25M_PCA';
%model_name = 'game';
%model_name = 'EMPA';
%model_name = 'VAE';
%model_name = 'VAE_e1k';
%what = 'conv3';
%what = 'linear2';
what = 'all';
%what = '';
%what = 'novelty';
%what = 'all';
%what = 'theory';
%what = 'termination';
%what = 'sprite';
project = 0;
glmodel = 1;
%suffix = '_nowhiten_nofilter';
%suffix = '_';
suffix = '_parts=123';
normalize = 1;
concat = 0;
novelty = 1;
saveYhat = 0;

%agg_filename = sprintf('mat/agg_gp_CV_us=%d_glm=%d_%s_fast.mat', use_smooth, glmodel, what);
%agg_filename = sprintf('mat/agg_gp_CV_us=%d_glm=%d_model=EMPA_%s_nsamples=100_fast_WTF.mat', use_smooth, glmodel, what);
%agg_filename = sprintf('mat/agg_gp_CV_us=%d_glm=%d_model=EMPA_%s_nsamples=100_fast.mat', use_smooth, glmodel, what);
%agg_filename = sprintf('/Volumes/fMRI-2/Mac_mat/agg_gp_CV_us=%d_glm=%d_model=%s_%s_nsamples=100_project=%d_fast=1.mat', use_smooth, glmodel, model_name, what, project);
% fasse
%agg_filename = fullfile(get_mat_dir(), sprintf('agg_gp_CV_us=%d_glm=%d_model=%s_%s_nsamples=100_project=%d_fast=1%s.mat', use_smooth, glmodel, model_name, what, project, suffix));
%agg_filename = fullfile(get_mat_dir(), sprintf('agg_gp_CV_us=%d_glm=%d_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=1%s.mat', use_smooth, glmodel, model_name, what, project, normalize, concat, novelty, suffix));
%agg_filename = fullfile(get_mat_dir(), sprintf('agg_gp_CV_cannon_repro_us=%d_glm=%d_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=1%s.mat', use_smooth, glmodel, model_name, what, project, normalize, concat, novelty, suffix));
agg_filename = fullfile(get_mat_dir(), sprintf('agg_gp_CV_us=%d_glm=%d_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=1%s.mat', use_smooth, glmodel, model_name, what, project, normalize, concat, novelty, suffix));
agg_filename

for s = 1:1:length(subjects)
    subj_id = subjects(s);

    %filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_%s_fast.mat', subj_id, use_smooth, glmodel, what);
    %filename = sprintf('mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=EMPA_%s_nsamples=100_fast.mat', subj_id, use_smooth, glmodel, what);
    %filename = sprintf('/Volumes/fMRI-2/Mac_mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=EMPA_%s_nsamples=100_fast.mat', subj_id, use_smooth, glmodel, what);
    %filename = sprintf('/Volumes/fMRI-2/Mac_mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=EMPA_%s_nsamples=100_project=%d_fast=1_.mat', subj_id, use_smooth, glmodel, what, project);
    %filename = sprintf('/Volumes/fMRI-2/Mac_mat/fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_nsamples=100_project=%d_fast=1.mat', subj_id, use_smooth, glmodel, model_name, what, project);
    %filename = fullfile(get_mat_dir(), sprintf('fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_nsamples=100_project=%d_fast=1.mat', subj_id, use_smooth, glmodel, model_name, what, project));
    %filename = fullfile(get_mat_dir(2), sprintf('fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=1_saveYhat=%d%s.mat', subj_id, use_smooth, glmodel, model_name, what, project, normalize, concat, novelty, saveYhat, suffix));
    %filename = fullfile(get_mat_dir(2), sprintf('fit_gp_CV_HRR_cannon_repro_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=1_saveYhat=%d%s.mat', subj_id, use_smooth, glmodel, model_name, what, project, normalize, concat, novelty, saveYhat, suffix));
    %filename = fullfile(get_mat_dir(2), sprintf('fit_gp_CV_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=1_saveYhat=%d%s.mat', subj_id, use_smooth, glmodel, model_name, what, project, normalize, concat, novelty, saveYhat, suffix)); % this should be it_gp_CV_HRR_cannon_repro_, for normalize=2
    %filename = fullfile(get_mat_dir(0), sprintf('fit_gp_CV_HRR_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_nsamples=100_project=%d_norm=%d_fast=1.mat', subj_id, use_smooth, glmodel, model_name, what, project, normalize)); % DQN normalize=2
    %filename = fullfile(get_mat_dir(2), sprintf('fit_gp_CV_subj=%d_us=%d_glm=%d_mask=mask_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=1_saveYhat=%d%s.mat', subj_id, use_smooth, glmodel, model_name, what, project, normalize, concat, novelty, saveYhat, suffix)); % DQN normalize=2, project=0
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=DQN25M_all_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0.mat', subj_id); % DQN 25M
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=DQN25M_all_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=0.mat', subj_id);
    %filename = sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=game_all_nsamples=100_project=0_norm=1_fast=1_saveYhat=0.mat', subj_id)
    %filename = sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=EMPA_theory_nsamples=100_project=0_norm=1_fast=1_saveYhat=0.mat', subj_id);
    %filename = sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=PCA_all_nsamples=100_project=0_norm=1_fast=1_saveYhat=0.mat', subj_id);
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=VAE__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0.mat', subj_id); % VAE e1k
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=VAE__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=0.mat', subj_id); % VAE e1k
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=1_parts=123.mat', subj_id);
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=12.mat', subj_id); % fit_gp_CV from paper, except partitions 1 and 2
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=23.mat', subj_id); % fit_gp_CV from paper, except partitions 2 and 3
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=game__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123.mat', subj_id); % game id ???????????? wtf..........
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=game__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=1_parts=123_games=all.mat', subj_id); % game id ???????????? wtf..........
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=state__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=1_parts=123_games=all.mat', subj_id); % state features, z-scored across time, and normalized across features
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=state__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=1_parts=123_games=all.mat', subj_id); % state features, z-scored across time, and normalized across features
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=DQN25M_PCA_all_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123_games=all.mat', subj_id); %  DQN + PCA
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=DQN25M_PCA_all_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123_games=all.mat', subj_id); %  DQN + PCA
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=DQN25M_PCA_conv3_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123_games=all.mat', subj_id); %  DQN + PCA
    %filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=DQN25M_PCA_linear2_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123_games=all.mat', subj_id); %  DQN + PCA
    filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=DQN25M_PCA_all_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_saveYhat=0_parts=123_games=all.mat', subj_id); %  DQN + PCA
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

