% predict bold signal using leave-one-run-out CV

close all;
clear all;

EXPT = vgdl_expt();

subjects = [1:14 16:length(EXPT.subject)];

[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');

load(fullfile(get_mat_dir(false), 'glm_bic_bms_atlas=AAL2_GLM_102_multiplex_with_controls.mat'), 'glms');
glmodels = glms;


for glmodel = glmodels

    fprintf('  glm %d\n', glmodel);

    for s = 1:length(subjects)
        subj_id = subjects(s);

        tic
        fprintf('           loading BOLD for %d...\n', subj_id);
        [Y, ~, ~, ~, run_id, ~, X, ~, KWX] = load_BOLD(EXPT, glmodel, subj_id, whole_brain_mask, Vwhole_brain_mask);
        toc

        tic
        fprintf('           loading betas for %d...\n', subj_id);
        [B, names] = ccnl_get_beta_series(EXPT, glmodel, subj_id, ' ', whole_brain_mask); % all
        B_CV = nan(size(B));

        for i = 1:length(names)
            suffix = names{i}(5:end);
            ix = find(contains(names, suffix) & ~strcmp(names, names{i}));
            B_CV(i,:) = mean(B(ix,:), 1); % average betas from other runs
        end

        %predict
        Yhat = X * B_CV;

        % Compute correlation separately for every run
        run_ids = unique(run_id);
        r = nan(length(run_ids), size(Yhat, 2));
        p = nan(length(run_ids), size(Yhat, 2));
        for k = 1:length(run_ids)
            which = (run_id == run_ids(k));
            for i = 1:size(Yhat, 2)
                [r(k,i), p(k,i)] = corr(Y(which,i), Yhat(which,i));
            end
        end
        
        filename = sprintf('glm_predict3_glm=%d_subj=%d.mat', glmodel, subj_id);
        fprintf('Saving to %s\n', filename);
        clear Y;
        clear Yhat;
        clear KWX;
        clear B;
        save(fullfile(get_mat_dir(2), filename));
    end
    
end
