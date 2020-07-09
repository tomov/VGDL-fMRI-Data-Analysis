% plot results from agg_gp_CV.m

clear;

load('mat/agg_gp_CV_us=1_theory.mat');


if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

%bspmview_wrapper(EXPT, tmap);
bspmview_wrapper(EXPT, logGBFmap);
