clear all;

load('mat/agg_neurosynth_rsa2_us=1_glm=9_model=EMPA_theory_project=0.mat');

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

%% odd subjects only
%subjs = 1:2:32;
%[h,p,ci,stats] = ttest(zs(subjs, :));
%%[h,p,ci,stats] = ttest(zs(subjs, :), zs_null(subjs, :));
%ts = stats.tstat;
%tmap(mask) = ts;


%bspmview_wrapper(EXPT, SMSE_tmap);
bspmview_wrapper(EXPT, tmap);
%bspmview_wrapper(EXPT, diff_R2map);
%bspmview_wrapper(EXPT, null_R2map);