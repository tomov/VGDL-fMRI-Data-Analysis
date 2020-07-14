% plot results from agg_gp_CV.m

close all;
clear all;

load('mat/agg_gp_CV_us=1_glm=9_theory.mat');
%load('mat/agg_gp_CV_us=1_sprite.mat');

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
