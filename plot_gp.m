% plot results from agg_gp.m

clear;

load('mat/agg_gp_us=1_theory.mat');

lme = -group_bic/2;

map = zeros(size(mask));
%map(mask) = mean_adjR2;
map(mask) = lme + max(abs(lme));
%map(mask) = mean_sigma;


if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

bspmview_wrapper(EXPT, map);
