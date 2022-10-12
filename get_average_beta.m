% get average beta from GLM across a bunch of subjects

close all;
clear all;

EXPT = vgdl_expt();


train_subjects = 1:2:length(EXPT.subject);

[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');

load(fullfile(get_mat_dir(false), 'glm_bic_bms_atlas=AAL2_GLM_102_multiplex_with_controls.mat'), 'glms');
glmodels = glms;

for glmodel = glmodels

    fprintf('  glm %d\n', glmodel);
    
    tic
    disp('           loading initial betas...');
    [B, names] = ccnl_get_beta_series(EXPT, glmodel, 1, ' ', whole_brain_mask); % all
    toc

    names_no_sess = [];
    for i = 1:length(names)
        names_no_sess{i} = get_unique_regressor_substring(names{i});
    end

    unique_regressors = unique(names_no_sess);
    cum_B = zeros(length(unique_regressors), size(B, 2));

    for s = 1:length(train_subjects)
        subj_id = train_subjects(s);

        tic
        fprintf('           loading betas for %d...\n', subj_id);
        [B, names] = ccnl_get_beta_series(EXPT, glmodel, subj_id, ' ', whole_brain_mask); % all
        toc

        for r = 1:length(unique_regressors)
            substring = unique_regressors{r};
            which = contains(names, substring);
            assert(sum(which) >= 1);
            mean_B = mean(B(which, :), 1);
            cum_B(r,:) = cum_B(r,:) + mean_B;
        end

    end

    avg_B = cum_B / length(train_subjects);

    filename = sprintf('get_average_beta_glm=%d.mat', glmodel);
    fprintf('Saving to %s\n', filename);
    save(fullfile(get_mat_dir(2), filename));
end
