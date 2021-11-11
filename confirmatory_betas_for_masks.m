% extract beta coefficients for masks
function confirmatory_betas_for_masks(glmodel, contrast, Num, sphere)

%clear all;

%glmodel = 102;
%contrast = 'theory_change_flag';
%Num = 1;
%sphere = 10;

%confirmatory_regressors = {'theory_change_flag', 'up', 'down', 'left', 'right', 'spacebar', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', 'avatar_collision_flag', 'effectsByCol', 'play_start', 'play_end'};
%confirmatory_glm = repmat(102, [1 length(confirmatory_regressors)]);
%
%confirmatory_regressors = [confirmatory_regressors, {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'}];
%confirmatory_glm = [confirmatory_glm 3 85 51 52];

confirmatory_regressors =  {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};
confirmatory_glm = [103 104 105];

EXPT = vgdl_expt();

[subj_ids, subjdirs, goodRuns, goodSubjects] = vgdl_getSubjectsDirsAndRuns();

% spherical mask around top ROI from contrast
if ischar(glmodel)
    if ismember(glmodel, {'AAL2', 'AAL3v1', 'HarvardOxford', 'AAL2_grouped3', 'AAL2_ungrouped', 'AAL2_grouped4', 'AAL2_ungrouped2'})
        % anatomical ROI
        atlas_name = glmodel;
        [mask_filenames, regions] = get_anatomical_masks(atlas_name);
        filename = fullfile(get_mat_dir(false), sprintf('confirmatory_betas_for_masks_atlas=%s_cglm=%s.mat', atlas_name, sprintf('%d-', confirmatory_glm)));
    else
        % a priori ROIs
        tag = glmodel; % fake "glmodel" = study tag
        [mask_filenames, regions] = get_masks_from_study(tag, sphere);
        filename = fullfile(get_mat_dir(false), sprintf('confirmatory_betas_for_masks_tag=%s_sphere=%.1fmm_cglm=%s.mat', tag, sphere, sprintf('%d-', confirmatory_glm)));
    end
    glmodel = 9; % for load_BOLD; doesn't really matter
else
    % actual GLM
    [mask_filenames, regions] = get_masks_from_contrast(glmodel, contrast, true, [], Num, sphere);
    filename = fullfile(get_mat_dir(false), sprintf('confirmatory_betas_for_masks_glm=%d_con=%s_Num=%d_sphere=%.1fmm_cglm=%s.mat', glmodel, contrast, Num, sphere, sprintf('%d-', confirmatory_glm)));
end
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

