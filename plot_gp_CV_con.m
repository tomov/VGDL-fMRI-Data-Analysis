% plot contrast results from agg_gp_CV.m

close all;
clear all;

%load('mat/agg_gp_CV_us=1_glm=9_theory.mat');
%load('mat/agg_gp_CV_us=1_sprite.mat');

load('mat/agg_gp_CV_us=1_sprite.mat');
zs_sprite = zs;

load('mat/agg_gp_CV_us=1_termination.mat');
zs_termination = zs;

load('mat/agg_gp_CV_us=1_interaction.mat');
zs_interaction = zs;

%[h,p,ci,stats] = ttest(zs_interaction, zs_sprite);
%[h,p,ci,stats] = ttest(zs_interaction, zs_termination);
%[h,p,ci,stats] = ttest(zs_termination, zs_sprite);
ts = stats.tstat;
tmap(mask) = ts;



if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end


bspmview_wrapper(EXPT, tmap);
