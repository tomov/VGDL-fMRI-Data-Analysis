clear all;

%load('mat/agg_ridge_CV_us=1_glm=9_theory.mat');
%load('mat/agg_ridge_CV_us=1_glm=9_subsample=1_theory.mat');
%load('mat/agg_ridge_CV_us=1_glm=9_model=EMPA_theory_subsample=0_project=1.mat');

%load('mat/agg_ridge_CV_us=1_glm=1_model=EMPA_theory_subsample=0_project=1.mat');

%load('mat/agg_ridge_CV_us=1_glm=9_model=game__subsample=0_project=0.mat');
%load('mat/agg_ridge_CV_us=1_glm=9_model=nuisance__subsample=0_project=0.mat');
load('mat/agg_ridge_CV_us=1_glm=101_model=nuisance__subsample=0_project=0.mat');
zs_null = zs;

load('mat/agg_ridge_CV_us=1_glm=9_model=EMPA_theory_subsample=0_project=0.mat');
%% odd subjects only
subjs = 1:2:32;
%[h,p,ci,stats] = ttest(zs(subjs, :));
[h,p,ci,stats] = ttest(zs(subjs, :), zs_null(subjs, :));
ts = stats.tstat;
tmap(mask) = ts;


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
