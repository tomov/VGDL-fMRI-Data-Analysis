% aggregate results from fit_gp_CV_noHRF.m
% copy of agg_gp_CV.m

use_smooth = true;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);
what = 'theory';
glmodel = 9;

%agg_filename = sprintf('mat/agg_gp_CV_noHRF_us=%d_glm=%d_%s_subsample=1.mat', use_smooth, glmodel, what);
agg_filename = sprintf('mat/agg_gp_CV_noHRF_us=%d_glm=%d_%s_subsample=0.mat', use_smooth, glmodel, what);

for s = 1:length(subjects)
    subj_id = subjects(s);

    %filename = sprintf('mat/fit_gp_CV_noHRF_HRR_subj=%d_us=%d_glm=%d_mask=mask_subsample=1_%s.mat', subj_id, use_smooth, glmodel, what);
    filename = sprintf('mat/fit_gp_CV_noHRF_HRR_subj=%d_us=%d_glm=%d_mask=mask_subsample=0_%s.mat', subj_id, use_smooth, glmodel, what);
    load(filename, 'sigmas', 'R2_CV', 'r_CV', 'MSE_CV', 'SMSE_CV', 'subj', 'use_smooth', 'glmodel', 'mask', 'what');

    % calc stuff
    k = 1; % 1 param = sigma
    n = 1698; % # TRs

    if s == 1
        rs = nan(length(subjects), size(r_CV,2));

        adjR2s = nan(length(subjects), size(r_CV,2));

        R2s = nan(length(subjects), size(r_CV,2));

        MSEs = nan(length(subjects), size(r_CV,2));
        SMSEs = nan(length(subjects), size(r_CV,2));
    end

    % MSEs
    MSEs(s,:) = mean(MSE_CV, 1);
    SMSEs(s,:) = mean(SMSE_CV, 1);
    % R2s
    R2s(s,:) = mean(R2_CV, 1);

    %adjR2s(s,:) = adjR2;

    % CV Pearson r's across subjects
    rs(s,:) = mean(r_CV, 1);

end

% Fisher z transform
zs = atanh(rs);

% t-test Pearson corr across subjects
[h,p,ci,stats] = ttest(zs);
%[h,p,ci,stats] = ttest(zs, null_zs);
ts = stats.tstat;



R2map = zeros(size(mask));
R2map(mask) = mean(R2s,1);

adjR2map = zeros(size(mask));
adjR2map(mask) = mean(adjR2s,1);

tmap = zeros(size(mask));
tmap(mask) = ts;


agg_filename
save(agg_filename);

