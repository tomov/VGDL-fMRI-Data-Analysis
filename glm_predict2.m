
close all;
clear all;

EXPT = vgdl_expt();

test_subjects = 2:2:length(EXPT.subject);

[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');

load(fullfile(get_mat_dir(false), 'glm_bic_bms_atlas=AAL2_GLM_102_multiplex_with_controls.mat'), 'glms');
glmodels = glms;


for glmodel = glmodels

    fprintf('  glm %d\n', glmodel);

    filename = sprintf('get_average_beta2_glm=%d.mat', glmodel);
    load(fullfile(get_mat_dir(2), filename), 'avg_names', 'avg_B');

    for s = 1:length(test_subjects)
        subj_id = test_subjects(s);

        tic
        fprintf('           loading BOLD for %d...\n', subj_id);
        [Y, ~, ~, ~, ~, ~, X, ~, KWX] = load_BOLD(EXPT, glmodel, subj_id, whole_brain_mask, Vwhole_brain_mask);
        toc

        tic
        fprintf('           loading betas for %d...\n', subj_id);
        [B, names] = ccnl_get_beta_series(EXPT, glmodel, subj_id, ' ', whole_brain_mask); % all
        B = nan(size(B));

        for i = 1:length(names)
            ix = find(strcmp(avg_names, names{i}));
            if isempty(ix)
                ix = find(contains(avg_names, names{i}(5:end)));
                B(i,:) = mean(avg_B(ix,:), 1);
            else
                B(i,:) = avg_B(ix,:);
            end
        end

        %predict
        Yhat = X * B;

        r = nan(1, size(Yhat, 2));
        p = nan(1, size(Yhat, 2));
        for i = 1:size(Yhat, 2)
            [r(i), p(i)] = corr(Y(:,i), Yhat(:,i));
        end
        
        filename = sprintf('glm_predict2_glm=%d_subj=%d.mat', glmodel, subj_id);
        fprintf('Saving to %s\n', filename);
        clear Y;
        clear Yhat;
        clear KWX;
        clear B;
        save(fullfile(get_mat_dir(2), filename));
    end
    
end
