% predict bold signal using leave-one-run-out CV, compute residuals

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
        fprintf('           loading SPM for %d...\n', subj_id);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj_id)]);
        load(fullfile(modeldir,'SPM.mat'));
        toc

        tic
        fprintf('           loading BOLD for %d...\n', subj_id);
        [~, ~, ~, ~, run_id, ~, X, KWY, KWX] = load_BOLD(EXPT, glmodel, subj_id, whole_brain_mask, Vwhole_brain_mask);
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
        %KWYhat = KWX * B; % see spm_spm.m
        KWYhat = KWX * B_CV; % see spm_spm.m

        % residuals
        res = KWY - KWYhat;
        ResSS = sum(res.^2,1);
        ResMS = ResSS / SPM.xX.trRV;

        %sanity check (use B not B_CV)
        %V = spm_vol(fullfile(modeldir,'ResMS.nii'));
        %ResMS1 = spm_data_read(V, whole_brain_mask);
        %immse(ResMS,ResMS1)
        
        filename = sprintf('glm_predict4_glm=%d_subj=%d.mat', glmodel, subj_id);
        fprintf('Saving to %s\n', filename);
        clear KWY;
        clear KWYhat;
        clear KWX;
        clear B;
        clear B_CV;
        clear res;
        save(fullfile(get_mat_dir(2), filename));
    end
    
end
