% RSA with ground truth theory HRRs, with patterns averaged across partitions
% do it manually b/c we average across partitions (runs), which is not supported by pipeline

%TODO doesn't work -- lots of negative corr's; probably due to within-run and across-run comparisons being mixed 

% copied from neurosynth_rsa_HRR and neurosynth_rsa_HRR_groundtruth


clear all;

use_smooth = false;
glmodel = 1;
nperms = 1000;
dist = 'euclidean';

radius = 6 / 1.5; % 5 mm

if use_smooth
    EXPT = vgdl_expt();
    [whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');
else
    EXPT = vgdl_expt_nosmooth();
    [whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask_nosmooth.nii');
end

subjects = 1:length(EXPT.subject);

% get neurosynth ROIs
% % [roi_masks, region] = get_neurosynth_rois(lateralized); too slow
load('mat/get_neurosynth_rois_lat=true'); % generated by get_neurosynth_rois
filename = sprintf('mat/neurosynth_rsa_HRR_groundtruth_avg_us=%d_nperms=%d_dist=%s.mat', use_smooth, nperms, dist);

% or, searchlight ROIs
% % ROI = get_searchlight_rois(whole_brain_mask, Vwhole_brain_mask, radius);  too slow
%load(sprintf('mat/get_searchlight_rois_us=%d_r=%.4f.mat', use_smooth, radius));
%roi_masks = {ROI.voxel_idx}; % backwards compatibility; they're indices actually, not masks, but it'll do
%filename = sprintf('mat/searchlight_rsa_HRR_groundtruth_avg_us=%d_r=%.4f_nperms=%d_dist=%s.mat', use_smooth, radius, nperms, dist);

%load('mat/HRR_groundtruth_RDM_correlation.mat'); % game_names, mean_RDM
load(sprintf('mat/HRR_groundtruth_RDM_K=10_N=10_E=0.050_nsamples=100_dist=%s.mat', dist)); % game_names, mean_RDM
game_names = cellfun(@strtrim, mat2cell(game_names, ones(size(game_names, 1), 1)), 'UniformOutput', false);
behavioral_RDM = mean_RDM;
ng = length(game_names);

% allocate memory
for r = 1:length(roi_masks)
    ROI(r).rho = nan(length(subjects), 1);
    ROI(r).null_rho = nan(length(subjects), nperms);
end

% for each subject
%
for s = 1:length(subjects)
    subj_id = subjects(s);
    subj_id

    disp('getting betas and RDMs...');
    tic

    % precompute whole-brain patterns (betas) for each game
    %
    for g = 1:ng
        game_name = game_names{g};

        B = ccnl_get_beta_series(EXPT, glmodel, subj_id, game_name, whole_brain_mask);
        assert(size(B,1) == 3);
        %B = B(:, ~any(isnan(B), 1));
        assert(~any(isnan(B(:))));

        % average across partitions (see Diedrichsen and Kriegeskorte, 2017)
        U_all(g,:) = mean(B, 1);
    end

    upper = logical(triu(ones(size(mean_RDM)), 1));

    toc
    
    disp('ROIs....');
    tic

    % for each ROI
    %
    for r = 1:length(roi_masks)
        roi_mask = roi_masks{r};

        if mod(r,1) == 0
            r
            toc
            tic
        end

        % subset whole-brain patterns
        U = U_all(:, roi_mask(whole_brain_mask));

        % compute neural RDM (ccnl_roi_rdms.m)
        neural_RDM = squareRDMs(pdist(U, dist));

        % match RDMs with rank correlation (ccnl_match_rdms.m)
        %
        assert(isequal(size(neural_RDM), size(mean_RDM)), 'Neural and behavioral RDM sizes should be equal');

        rho = corr(neural_RDM(upper), mean_RDM(upper), 'type', 'Spearman');

        ROI(r).rho(s) = rho;

        % null distr
        %
        for i = 1:nperms
            U = U(randperm(ng), :);

            neural_RDM = squareRDMs(pdist(U, dist));
            null_rho = corr(neural_RDM(upper), mean_RDM(upper), 'type', 'Spearman');

            ROI(r).null_rho(s,i) = null_rho;
        end
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

    % note the importance of comparing aggregate statistics to aggregate statistics
    ROI(r).avg_rho = mean(ROI(r).rho, 1); % across subj
    ROI(r).avg_null_rho = [mean(ROI(r).null_rho, 1) ROI(r).avg_rho]; % include actual rho in null
    ROI(r).null_p = mean(ROI(r).avg_rho < ROI(r).avg_null_rho);
    ROI(r).null_z = (ROI(r).avg_rho - mean(ROI(r).avg_null_rho)) / std(ROI(r).avg_null_rho);
end

save(filename, '-v7.3');

