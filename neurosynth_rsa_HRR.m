% RSA with ground truth theory HRRs
% do it manually b/c we average across partitions (runs), which is not supported by pipeline

use_smooth = false;
glmodel = 1;


if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);


% get ROIs
%[roi_masks, region] = get_neurosynth_rois(lateralized);
load('mat/get_neurosynth_rois_lat=true');
mask = roi_masks{region == 120};

load('mat/HRR_groundtruth_RDM_correlation.mat'); % game_names, mean_RDM
game_names = cellfun(@strtrim, mat2cell(game_names, ones(size(game_names, 1), 1)), 'UniformOutput', false);

behavioral_RDM = mean_RDM;

upper = logical(triu(ones(size(neural_RDM)), 1));
lower = logical(tril(ones(size(neural_RDM)), -1));

clear S;

for s = 1:length(subjects)
    subj = subjects(s);

    clear U;
    for g = 1:length(game_names)
        game_name = game_names{g};

        B = ccnl_get_beta_series(EXPT, glmodel, subj, game_name, mask);
        assert(size(B,1) == 3);
        B = B(:, ~any(isnan(B), 1));

        U(g,:) = B(2,:);
        U(g + length(game_names),:) = B(3,:);
    end

    % RDM between partition 2 and partition 3
    % 1) avoid comparing games in the same run
    % 2) look at more asymptotic behavior (partition 1 is probs mostly learning)
    tmp = squareRDMs(pdist(U, 'correlation'));
    neural_RDM = tmp(1:length(game_names), length(game_names)+1:end);

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
    S(s).rho = corr(neural_RDM(upper), behavioral_RDM(upper), 'type', 'Spearman');
    S(s).z_rho = atanh(S(s).rho); % should be high

end


sym = [S.sym];
diff = [S.diff];
z_rho = [S.z_rho];

figure;

%sem = @(x) std(x, 0, 1) / sqrt(size(x,1));
sem = @(x) std(x) / sqrt(length(x));

m = [mean(sym) mean(diff) mean(z_rho)];
se = [sem(sym) sem(diff) sem(z_rho)];
bar([1 2 3], m);
hold on;
er = errorbar([1 2 3], m, se, 'linestyle', 'none');

