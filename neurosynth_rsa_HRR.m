% RSA with ground truth theory HRRs
% do it manually b/c we average across partitions (runs), which is not supported by pipeline

use_smooth = true;
glmodel = 1;
nperms = 1000;

filename = sprintf('mat/neurosynth_rsa_HRR_us=%d_nperm=%d.mat', use_smooth, nperms);

if use_smooth
    EXPT = vgdl_expt();
    whole_brain_mask = ccnl_load_mask('masks/mask.nii');
else
    EXPT = vgdl_expt_nosmooth();
    whole_brain_mask = ccnl_load_mask('masks/mask_nosmooth.nii');
end

subjects = 1:length(EXPT.subject);


% get ROIs
%[roi_masks, region] = get_neurosynth_rois(lateralized);
load('mat/get_neurosynth_rois_lat=true');

load('mat/HRR_groundtruth_RDM_correlation.mat'); % game_names, mean_RDM
game_names = cellfun(@strtrim, mat2cell(game_names, ones(size(game_names, 1), 1)), 'UniformOutput', false);

behavioral_RDM = mean_RDM;

upper = logical(triu(ones(size(behavioral_RDM)), 1));
lower = logical(tril(ones(size(behavioral_RDM)), -1));

%mask = roi_masks{region == 120};

% precompute whole-brain patterns for each subject
%
clear U_all;

for s = 1:length(subjects)
    subj = subjects(s);
    for g = 1:length(game_names)
        game_name = game_names{g};

        B = ccnl_get_beta_series(EXPT, glmodel, subj, game_name, whole_brain_mask);
        assert(size(B,1) == 3);
        %B = B(:, ~any(isnan(B), 1));
        assert(~any(isnan(B(:))));

        %{
        % sanity
        B1 = ccnl_get_beta_series(EXPT, glmodel, subj, game_name, mask);
        B1 = B1(:, ~any(isnan(B1), 1));
        B2 = B(:,mask(whole_brain_mask));
        assert(immse(B1, B2) < 1e-10);
        %}

        % for each game, compare patterns from partition 2 vs. 3
        U_all{s}(g,:) = B(2,:);
        U_all{s}(g + length(game_names),:) = B(3,:);
    end
end


clear ROI;

% for each ROI
%
for r = 1:length(roi_masks)
    roi_mask = roi_masks{r};

    clear S;

    tic

    % 
    for s = 1:length(subjects)
        subj = subjects(s);

        % subset whole-brain patterns
        U = U_all{s}(:, whole_brain_mask(roi_mask));

        % RDM between partition 2 and partition 3
        % 1) avoid comparing games in the same run
        % 2) look at more asymptotic behavior (partition 1 is probs mostly learning)
        neural_RDM = 1 - corr(U(1:length(game_names),:)', U(length(game_names)+1:end,:)'); %-- alternative; not much faster

        %{
        % sanity -- note corr inputs are transposed compared to pdist 
        tmp = squareRDMs(pdist(U, 'correlation'));
        neural_RDM1 = tmp(1:length(game_names), length(game_names)+1:end);
        assert(immse(neural_RDM1, neural_RDM) < 1e-10);
        %}

        % how symmetric is the RDM?
        S(s).sym = how_sym(neural_RDM); % should be high

        S(s).neural_RDM = neural_RDM;

        % make it symmetric by averaging two halves
        neural_RDM = 0.5 * (neural_RDM + neural_RDM'); 

        % similarity of same vs different games
        S(s).within_game = mean(diag(neural_RDM));
        S(s).across_game = mean(neural_RDM(upper));
        S(s).diff = S(s).across_game - S(s).within_game; % should be high

        % rank correlation between behavioral and neural RDMs
        S(s).rho = corr(neural_RDM(upper), behavioral_RDM(upper), 'type', 'Spearman'); % should be high

        S(s).z_rho = atanh(S(s).rho);

        % null distribution TODO dedupe
        %
        S(s).null_sym = nan(nperms, 1);
        S(s).null_diff = nan(nperms, 1);
        S(s).null_rho = nan(nperms, 1);

        for i = 1:nperms
            U(1:length(game_names),:) = U(randperm(length(game_names)), :); % shuffle one half only

            neural_RDM = 1 - corr(U(1:length(game_names),:)', U(length(game_names)+1:end,:)');

            S(s).null_sym(i) = how_sym(neural_RDM); % should be high

            neural_RDM = 0.5 * (neural_RDM + neural_RDM'); 

            within_game = mean(diag(neural_RDM));
            across_game = mean(neural_RDM(upper));
            S(s).null_diff(i) = across_game - within_game; % should be high

            S(s).null_rho(i) = corr(neural_RDM(upper), behavioral_RDM(upper), 'type', 'Spearman');
        end

    end

    ROI(r).S = S;

    % avg across subjects
    sym = mean([S.sym]);
    diff = mean([S.diff]);
    rho = mean([S.rho]);

    % avg across subjects
    null_sym = [mean([S.null_sym], 2); sym];
    null_diff = [mean([S.null_diff], 2); diff];
    null_rho = [mean([S.null_rho], 2); rho];

    % note the importance of comparing aggregate statistics to aggregate statistics
    % it would be incorrect e.g. to generate a null of individual subject symmetries 
    % and then compare the average symmetry to it
    ROI(r).p_sym = mean(sym < null_sym);
    ROI(r).p_diff = mean(diff < null_diff);
    ROI(r).p_rho = mean(rho < null_rho);

    [~, ROI(r).p_sym_t] = ttest([S.sym], 0, 'tail', 'right');
    [~, ROI(r).p_diff_t] = ttest([S.diff], 0, 'tail', 'right');
    [~, ROI(r).p_rho_t] = ttest([S.rho], 0, 'tail', 'right');

    ROI(r).z_sym = (sym - mean(null_sym)) / std(null_sym);
    ROI(r).z_diff = (diff - mean(null_diff)) / std(null_diff);
    ROI(r).z_rho = (rho - mean(null_rho)) / std(null_rho);

    ROI(r).sym = sym;
    ROI(r).diff = diff;
    ROI(r).rho = rho;

    ROI(r).null_sym = null_sym;
    ROI(r).null_diff = null_diff;
    ROI(r).null_rho = null_rho;

    toc
end

save(filename, '-v7.3');
