% plot results from agg_gp_CV.m

close all;
clear all;

%load('mat/agg_gp_CV_us=1_glm=9_model=EMPA_theory_nsamples=1_fast.mat');

%load('mat/agg_gp_CV_us=1_glm=9_model=EMPA_theory_nsamples=100_fast.mat');
%load('mat/agg_gp_CV_us=1_glm=9_model=EMPA_termination_nsamples=100_fast.mat');
%load('mat/agg_gp_CV_us=1_glm=9_model=EMPA_sprite_nsamples=100_fast.mat');
%load('mat/agg_gp_CV_us=1_glm=9_model=EMPA_interaction_nsamples=100_fast.mat');

%load('mat/agg_gp_CV_us=1_glm=9_theory_nsamples=100_project=0_fast=1_odd.mat');
%load('mat/agg_gp_CV_us=1_glm=9_theory_nsamples=100_project=0_fast=1_even.mat');
%load('mat/agg_gp_CV_us=1_glm=9_theory_nsamples=100_project=0_fast=1.mat');
%load('mat/agg_gp_CV_us=1_glm=9_sprite_nsamples=100_project=0_fast=1.mat');
%load('mat/agg_gp_CV_us=1_glm=9_interaction_nsamples=100_project=0_fast=1.mat');
%load('mat/agg_gp_CV_us=1_glm=9_termination_nsamples=100_project=0_fast=1.mat');

%load('mat/agg_gp_CV_us=1_glm=1_theory_nsamples=100_project=1_fast=1.mat');
%load('mat/agg_gp_CV_us=1_glm=1_model=nuisance__nsamples=100_project=0_fast=1.mat');

%% odd subject only
%subjs = 1:2:32;
%[h,p,ci,stats] = ttest(zs(subjs,:));
%ts = stats.tstat;
%tmap(mask) = ts;


%% t-test Pearson corr across subjects
%
%load('mat/agg_gp_CV_us=1_glm=1_model=game__nsamples=100_project=0_fast=1.mat');
%load('mat/agg_gp_CV_us=1_glm=1_model=nuisance__nsamples=100_project=0_fast=1.mat');

% paired t-test
%{
%load('mat/agg_gp_CV_us=1_glm=101_model=nuisance__nsamples=100_project=0_fast=1.mat');
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_all_nsamples=100_project=1_norm=1_fast=1.mat')); % !!!
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=VAE__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1.mat'));
load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=DQN25M_all_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1.mat');
zs_control = zs;
%load('mat/agg_gp_CV_us=1_glm=9_theory_nsamples=100_project=0_fast=1.mat');
load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1.mat')); % this is it !!!!!!!!!!!!!!!!!!!!ontroller five
[h,p,ci,stats] = ttest(zs, zs_control);
ts = stats.tstat;
tmap(mask) = ts;
%}


%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_conv1_nsamples=100_project=1_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_conv2_nsamples=100_project=1_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_conv3_nsamples=100_project=1_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_linear1_nsamples=100_project=1_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_linear2_nsamples=100_project=1_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_all_nsamples=100_project=1_norm=1_fast=1.mat')); % !!!
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=DQN25M_all_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=DQN25M_all_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA_25M_e1k.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=PCA_all_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1.mat');


%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=PCA_all_nsamples=100_project=1_norm=1_fast=1.mat')); % !!!

%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_interaction_nsamples=100_project=1_fast=1_nowhiten_nofilter.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1_.mat')); % this is it
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1__odd.mat')); % t
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1__even.mat')); % t
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=9_model=EMPA_theory_nsamples=100_project=1_fast=1_nowhiten_nofilter.mat'));

%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1.mat')); % this is it !!!!!!!!!!!!!!!!!!!!ontroller five
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_sprite_nsamples=100_project=1_fast=1.mat')); % this is it
%$load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_interaction_nsamples=100_project=1_fast=1.mat')); % this is it
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_termination_nsamples=100_project=1_fast=1.mat')); % this is it


%load(fullfile(get_mat_dir(), 'agg_gp_CV_cannon_repro_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1.mat'));

% normalize=2
%load(fullfile(get_mat_dir(), 'agg_gp_CV_cannon_repro_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_norm=2_concat=0_novelty=1_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=DQN_all_nsamples=100_project=1_norm=2_concat=0_novelty=1_fast=1.mat'));
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=0_norm=2_concat=0_novelty=1_fast=1.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=DQN_all_nsamples=100_project=0_norm=2_concat=0_novelty=1_fast=1.mat');

%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=state__nsamples=100_project=1_norm=1_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=state__nsamples=100_project=0_norm=1_fast=1.mat'));  % !!!
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=irrelevant__nsamples=100_project=1_norm=1_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=irrelevant__nsamples=100_project=0_norm=1_fast=1.mat'));

% novelty
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=0_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_sprite_nsamples=100_project=1_norm=1_concat=0_novelty=0_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_interaction_nsamples=100_project=1_norm=1_concat=0_novelty=0_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_termination_nsamples=100_project=1_norm=1_concat=0_novelty=0_fast=1.mat'));
%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_novelty_nsamples=100_project=1_norm=1_concat=0_novelty=0_fast=1.mat'));

%load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=VAE__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1.mat'));
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=game__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=state__nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=state__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=123_by_game.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=DQN25M_PCA_all_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_us=1_glm=1_model=DQN25M_PCA_linear2_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');

%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_by_game_us=1_glm=1_model=EMPA_theory_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_parts=123_by_game.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_by_game_us=1_glm=1_model=DQN25M_PCA_all_nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_by_game_us=1_glm=1_model=PCA__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');
%load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_by_game_us=1_glm=1_model=PCA__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');
load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/agg_gp_CV_by_game_us=1_glm=1_model=VAE__nsamples=100_project=0_norm=1_concat=0_novelty=1_fast=1_parts=123.mat');

%load('mat/agg_gp_CV_us=1_glm=9_theory_fast.mat');
%load('mat/agg_gp_CV_us=1_glm=9_theory_nsamples=1_fast.mat');
%load('mat/agg_gp_CV_noHRF_us=1_glm=9_theory_subsample=1.mat');
%load('mat/agg_gp_CV_noHRF_us=1_glm=9_theory_subsample=0.mat');

%load('mat/agg_gp_CV_us=1_glm=9_theory.mat');
%load('mat/agg_gp_CV_us=1_glm=9_sprite.mat');
%load('mat/agg_gp_CV_us=1_glm=9_interaction.mat');
%load('mat/agg_gp_CV_us=1_glm=9_termination.mat');


%load('mat/agg_gp_CV_noRKW_us=1_glm=21_theory.mat');
%load('mat/agg_gp_CV_us=1_sprite.mat');
%load('mat/agg_gp_CV_us=1_glm=21_theory.mat');
%load('mat/agg_gp_CV_us=1_glm=21_sprite.mat');
%load('mat/agg_gp_CV_us=1_glm=21_interaction.mat');
%load('mat/agg_gp_CV_us=1_glm=21_termination.mat');

%{
load('mat/agg_gp_CV_us=1_termination.mat');
logGBFmap_termination = logGBFmap;

load('mat/agg_gp_CV_us=1_interaction.mat');
logGBFmap_interaction = logGBFmap;

logGBFmap = logGBFmap_interaction - logGBFmap_termination;
%}

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end


%bspmview_wrapper(EXPT, SMSE_tmap);
bspmview_wrapper(EXPT, tmap);
%bspmview_wrapper(EXPT, diff_R2map);
%bspmview_wrapper(EXPT, null_R2map);
%bspmview_wrapper(EXPT, R2map);
%bspmview_wrapper(EXPT, adjR2map);
%bspmview_wrapper(EXPT, Wmap);
%bspmview_wrapper(EXPT, logGBFmap);
%bspmview_wrapper(EXPT, loglikmap);
%bspmview_wrapper(EXPT, null_loglikmap);
%bspmview_wrapper(EXPT, null_loglikmap);
