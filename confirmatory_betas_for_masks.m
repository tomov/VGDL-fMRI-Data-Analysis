% extract beta coefficients for masks

clear all;

glmodel = 21;
contrast = 'theory_change_flag';
Num = 1;
sphere = 10;

%confirmatory_glm = [3 85 51 52];
%confirmatory_regressors = {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

confirmatory_regressors = {'theory_change_flag', 'up', 'down', 'left', 'right', 'spacebar', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', 'avatar_collision_flag', 'effectsByCol', 'block_start', 'block_end', 'instance_start', 'instance_end', 'play_start', 'play_end'};

confirmatory_glm = repmat(21, [1 length(confirmatory_regressors)]);

EXPT = vgdl_expt();

[subj_ids, subjdirs, goodRuns, goodSubjects] = vgdl_getSubjectsDirsAndRuns();

% spherical mask around top ROI from contrast
[mask_filenames, regions] = get_masks_from_contrast(glmodel, contrast, true, [], Num, sphere);
filename = sprintf('mat/confirmatory_betas_for_masks_glm=%d_con=%s_Num=%d_sphere=%.1fmm_cglm=%s.mat', glmodel, contrast, Num, sphere, sprintf('%d-', confirmatory_glm));
disp(filename);

for m = 1:length(mask_filenames)
    mask_filename = mask_filenames{m};
    [~, mask_name{m}, ~] = fileparts(mask_filename);
    disp(mask_name{m});

    for s = 1:length(subj_ids)
        subj_id = subj_ids(s);
        fprintf('Subject %d\n', subj_id);

        for i = 1:length(confirmatory_glm)
            B = ccnl_get_beta_series(EXPT, confirmatory_glm(i), subj_id, confirmatory_regressors{i}, mask_filename);
            %B = rand(sum(goodRuns{subj_id}), 124); % for debugging

            betas{m}(s,i) = mean(B(:));
        end
    end
end

filename
save(filename);

