% RSA with theory HRRs
% do it manually b/c we average RDMs across features, which is not supported by pipeline

clear all;

use_smooth = false;
glmodel = 24;
rsa_idx = 5;
%nperms = 1000;

filename = sprintf('mat/neurosynth_rsa_HRR_us=%d.mat', use_smooth);

if use_smooth
    EXPT = vgdl_expt();
    whole_brain_mask = ccnl_load_mask('masks/mask.nii');
else
    EXPT = vgdl_expt_nosmooth();
    whole_brain_mask = ccnl_load_mask('masks/mask_nosmooth.nii');
end

subjects = 1:length(EXPT.subject);

rsa = EXPT.create_rsa;
assert(rsa.glmodel == glmodel);


% get ROIs
%[roi_masks, region] = get_neurosynth_rois(lateralized);
load('mat/get_neurosynth_rois_lat=true');

%load('mat/HRR_groundtruth_RDM_correlation.mat'); % game_names, mean_RDM
game_names = cellfun(@strtrim, mat2cell(game_names, ones(size(game_names, 1), 1)), 'UniformOutput', false);

for s = 1:length(subjects)
    subj = subjects(s);

end


Behavioral = ccnl_behavioral_rdms(EXPT, rsa_idx, subjects);
assert(isequal(size(Behavioral(1).subj(1).RDM, subj_theory_RDMs{1})));

clear ROI;

% for each subject
%
for s = 1:length(subjects)
    subj = subjects(s);
    subj

    disp('getting betas...');
    tic

    % get HRR RDM (in lieu of ccnl_behavioral_rdms.m)
    %
    glmodel_orig = glmodel;
    load(sprintf('mat/HRR_subject_RDM_subj=%d_K=10_N=10_E=0.050_nsamples=10_dist=corr.mat', subj), 'theory_RDM');
    assert(glmodel == glmodel_orig);
    assert(dist == 'correlation');

    % precompute whole-brain patterns (betas) for each subject
    %
    B = ccnl_get_beta_series(EXPT, glmodel, subj, 'run_', whole_brain_mask);
    assert(size(B,1) == size(theory_RDM, 1);
    assert(~any(isnan(B(:)))); % b/c of mask

    % subset RDMs
    upper = logical(triu(ones(size(neural_RDM)), 1));
    upper = upper & Behavioral(1).subj(s).partition_RDM; % don't compare within the same partition (run)! BOLD is VERY autocorrelated
    theory_RDM = theory_RDM(upper);

    toc
    
    disp('ROIs....');
    tic

    % for each ROI
    %
    for r = 1:length(roi_masks)
        roi_mask = roi_masks{r};

        % subset whole-brain patterns
        U = B(:, whole_brain_mask(roi_mask));

        % compute neural RDM (ccnl_roi_rdms.m)
        neural_RDM = squareRDMs(pdist(U, 'correlation'));

        % match RDMs with rank correlation (ccnl_match_rdms.m)
        %
        assert(isequal(size(neural_RDM), size(theory_RDM)), 'Neural and behavioral RDMs should be equal -- check rsa.which_trials');

        rho = corr(neural_RDM, theory_RDM, 'type', 'Spearman');

        ROI(r).rho(s) = rho;
    end

    toc

end

save(filename, '-v7.3');

% group-level stats
%
for r = 1:length(roi_masks)

    ROI(r).F = atanh(ROI(r).rho);
    [h,p,ci,stat] = ttest(ROI(r).F);
    ROI(r).t = stat.tstat;
    ROI(r).p = p;
end

save(filename, '-v7.3');


