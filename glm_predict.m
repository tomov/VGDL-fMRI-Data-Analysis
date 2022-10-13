
close all;
clear all;

EXPT = vgdl_expt();

test_subjects = 2:2:length(EXPT.subject);

[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');

load(fullfile(get_mat_dir(false), 'glm_bic_bms_atlas=AAL2_GLM_102_multiplex_with_controls.mat'), 'glms');
glmodels = glms;


for glmodel = glmodels

    fprintf('  glm %d\n', glmodel);

    filename = sprintf('get_average_beta_glm=%d.mat', glmodel);
    load(fullfile(get_mat_dir(2), filename), 'unique_regressors', 'avg_B');

    for s = 1:length(test_subjects)
        subj_id = test_subjects(s);

        tic
        fprintf('           loading BOLD for %d...\n', subj_id);
        [~, ~, ~, ~, ~, ~, X, KWY, KWX] = load_BOLD(EXPT, glmodel, subj_id, whole_brain_mask, Vwhole_brain_mask);
        toc

        tic
        fprintf('           loading betas for %d...\n', subj_id);
        [B, names] = ccnl_get_beta_series(EXPT, glmodel, subj_id, ' ', whole_brain_mask); % all
        B = nan(size(B));

        for i = 1:length(names)
            ur = find(strcmp(unique_regressors, get_unique_regressor_substring(names{i})));
            B(i,:) = avg_B(ur,:);
        end

        %predict
        KWYhat = KWX * B;

        r = nan(1, size(KWYhat, 2));
        p = nan(1, size(KWYhat, 2));
        for i = 1:size(KWYhat, 2)
            [r(i), p(i)] = corr(KWY(:,i), KWYhat(:,i));
        end
        
        filename = sprintf('glm_predict_glm=%d_subj=%d_KWY.mat', glmodel, subj_id);
        fprintf('Saving to %s\n', filename);
        %clear Y;
        %clear Yhat1;
        clear KWYhat;
        clear KWX;
        clear KWY;
        clear B;
        save(fullfile(get_mat_dir(2), filename));
    end
    
end
