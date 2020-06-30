% RSA with theory HRRs
% do it manually b/c we average RDMs across features, which is not supported by pipeline

rmpath('/n/sw/helmod/apps/centos7/Core/spm/12.7487-fasrc01/external/fieldtrip/external/stats/'); % tcdf

clear all;

use_smooth = false;
glmodel = 24;
rsa_idx = 5;
%nperms = 1000;

radius = 6 / 1.5; % mm -> voxels

if use_smooth
    EXPT = vgdl_expt();
    [whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');
else
    EXPT = vgdl_expt_nosmooth();
    [whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask_nosmooth.nii');
end

subjects = 1:length(EXPT.subject);

what = 'theory';

%subjects = [1, 2];

% get neurosynth ROIs
% % [roi_masks, region] = get_neurosynth_rois(lateralized); too slow
%load('mat/get_neurosynth_rois_lat=true'); % generated by get_neurosynth_rois
%filename = sprintf('mat/neurosynth_rsa_HRR_us=%d_%s.mat', use_smooth, what);

% or, searchlight ROIs
% % ROI = get_searchlight_rois(whole_brain_mask, Vwhole_brain_mask, radius);  too slow
load(sprintf('mat/get_searchlight_rois_us=%d_r=%.4f.mat', use_smooth, radius));
roi_masks = {ROI.voxel_idx}; % backwards compatibility; they're indices actually, not masks, but it'll do
filename = sprintf('mat/searchlight_rsa_HRR_us=%d_r=%.4f_%s.mat', use_smooth, radius, what);

%load('mat/HRR_groundtruth_RDM_correlation.mat'); % game_names, mean_RDM
%game_names = cellfun(@strtrim, mat2cell(game_names, ones(size(game_names, 1), 1)), 'UniformOutput', false);

Behavioral = ccnl_behavioral_rdms(EXPT, rsa_idx, subjects);

% for each subject
%
for s = 1:length(subjects)
    subj_id = subjects(s);
    subj_id

    rsa = EXPT.create_rsa(rsa_idx, subj_id);
    assert(rsa.glmodel == glmodel);

    disp('getting betas and RDMs...');
    tic

    % get HRR RDM (in lieu of ccnl_behavioral_rdms.m)
    % generated by HRR.py in py_vgdl
    %
    glmodel_orig = glmodel;
    [~, hostname] = system('hostname');
    if ~isempty(strfind(hostname, 'omchil')) || ~isempty(strfind(hostname, 'dhcp-'))
        % local
        load(sprintf('mat/HRR_subject_RDM_subj=%d_K=10_N=10_E=0.050_nsamples=100_dist=correlation.mat', subj_id), 'theory_RDM', 'sprite_RDM', 'interaction_RDM', 'termination_RDM', 'dist', 'glmodel'); % local
    else
        % ncf
        load(sprintf('../../VGDL/py_vgdl/mat/HRR_subject_RDM_subj=%d_K=10_N=10_E=0.050_nsamples=100_dist=correlation.mat', subj_id), 'theory_RDM', 'sprite_RDM', 'interaction_RDM', 'termination_RDM', 'dist', 'glmodel');
    end
    assert(glmodel == glmodel_orig);
    assert(isequal(dist, 'correlation'));
    assert(isequal(size(Behavioral(1).subj(s).RDM), size(theory_RDM)));
    assert(isequal(size(theory_RDM), size(sprite_RDM)));
    assert(isequal(size(theory_RDM), size(interaction_RDM)));
    assert(isequal(size(theory_RDM), size(termination_RDM)));

    % precompute whole-brain patterns (betas) for each subject
    %
    B = ccnl_get_beta_series(EXPT, glmodel, subj_id, 'run_', whole_brain_mask);
    assert(size(B,1) == size(theory_RDM, 1));
    assert(~any(isnan(B(:)))); % b/c of mask

    % subset RDMs
    upper = logical(triu(ones(size(theory_RDM)), 1));
    upper = upper & Behavioral(1).subj(s).partition_RDM; % don't compare within the same partition (run)! BOLD is VERY autocorrelated

    switch what
        case 'theory'
            HRR_RDM = theory_RDM;
        case 'sprite'
            HRR_RDM = sprite_RDM;
        case 'interaction'
            HRR_RDM = interaction_RDM;
        case 'termination'
            HRR_RDM = termination_RDM;
        otherwise
            assert(false);
    end

    toc
    
    disp('ROIs....');
    tic

    % for each ROI
    %
    for r = 1:length(roi_masks)
        roi_mask = roi_masks{r};

        if mod(r,1000) == 0
            r
            toc
            tic
        end

        % subset whole-brain patterns
        U = B(:, roi_mask(whole_brain_mask));

        % compute neural RDM (ccnl_roi_rdms.m)
        neural_RDM = squareRDMs(pdist(U, 'correlation'));

        % match RDMs with rank correlation (ccnl_match_rdms.m)
        %
        assert(isequal(size(neural_RDM), size(HRR_RDM)), 'Neural and behavioral RDM sizes should be equal -- check rsa.which_trials');

        rho = corr(neural_RDM(upper), HRR_RDM(upper), 'type', 'Spearman');
        %[taub, ~, ~, ~, tauz] = ktaub([neural_RDM(upper) HRR_RDM(upper)], 0.05); too slow; untested

        ROI(r).rho(s) = rho;
        %ROI(r).taub(s) = taub;
        %ROI(r).tauz(s) = tauz;
    end

    toc

end

save(filename, '-v7.3');

% group-level stats
%
for r = 1:length(roi_masks)

    % t-test Spearman rho's
    ROI(r).F = atanh(ROI(r).rho);
    [h,p,ci,stat] = ttest(ROI(r).F);
    ROI(r).t_s = stat.tstat;
    ROI(r).p_s = p;

    % t-test Kendall tau's (standardized)
    %[h,p,ci,stat] = ttest(ROI(r).tauz);
    %ROI(r).t_k = stat.tstat;
    %ROI(r).p_k = p;
end

save(filename, '-v7.3');


