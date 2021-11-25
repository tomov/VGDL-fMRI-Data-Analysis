function gen_PETHs(glmodel, contrast, Num, sphere, what)
    % get activations around given event(s)

    EXPT = vgdl_expt();

    %subj_ids = [1:32];
    [subj_ids, subjdirs, goodRuns, goodSubjects] = vgdl_getSubjectsDirsAndRuns();

    %sphere = 4; % sphere radius in mm
    %Num = 3; % # peaks per ROI
    %glmodel = 21;
    %contrast = 'theory_change_flag';

    [~, whole_brain_mask, ~] = get_mask_format_helper('masks/mask.nii');

    PETH_dTRs = -2:10; % TRs relative to the event onset to use for the PETH's
    baseline_dTRs = -2:0; % TRs relative to the event onset to use for the baseline

    % spherical mask around top ROI from contrast

    if ischar(glmodel)
        if ismember(glmodel, {'AAL2', 'AAL3v1', 'HarvardOxford', 'AAL2_GP_EMPA_grouped', 'AAL2_GP_EMPA', 'AAL2_GLM_102_grouped', 'AAL2_GLM_102', 'AAL2_GP_EMPA_GLM_102_grouped', 'AAL2_GP_EMPA_GLM_102'})
            % anatomical ROI
            atlas_name = glmodel;
            [mask_filenames, regions] = get_anatomical_masks(atlas_name);
            filename = fullfile(get_mat_dir(false), sprintf('PETHs_atlas=%s_%s.mat', atlas_name, what));
        else
            % a priori ROIs
            tag = glmodel; % fake "glmodel" = study tag
            [mask_filenames, regions] = get_masks_from_study(tag, sphere);
            filename = fullfile(get_mat_dir(false), sprintf('PETHs_tag=%s_sphere=%.1fmm_%s.mat', tag, sphere, what));
        end
        glmodel = 9; % for load_BOLD; doesn't really matter
    else
        % actual GLM
        [mask_filenames, regions] = get_masks_from_contrast(glmodel, contrast, true, [], Num, sphere);
        filename = fullfile(get_mat_dir(false), sprintf('PETHs_glm=%d_con=%s_odd_Num=%d_sphere=%.1fmm_%s.mat', glmodel, contrast, Num, sphere, what));
    end
    disp(filename);

    % which events to extract time courses for
    regs_fields = {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};
    visuals_fields = {'effects', 'avatar_collision_flag', 'new_sprites', 'killed_sprites'};
    onoff_fields = {'block_start', 'block_end', 'instance_start', 'instance_end', 'play_start', 'play_end'};
    fields = [regs_fields, visuals_fields, onoff_fields];

    % where we store the timecourses, for averaging in the end
    activations = struct;
    counts = struct;

    % loop over subjects
    for s = 1:length(subj_ids)
        subj_id = subj_ids(s);

        fprintf('Subject %d\n', subj_id);

        % event onsets, by type
        onsets = struct;

        disp('extracting event onsets');
        tic

        % extract event onsets for each run
        for SPM_run_id = 1:length(EXPT.subject(subj_id).functional)

            run_id = get_behavioral_run_id(subj_id, SPM_run_id);

            % in lieu of get_regressors, get_visuals, get_onoffs, etc.
            % because there's no mongo on NCF
            load(fullfile(get_mat_dir(true), sprintf('get_regressors_subj%d_run%d.mat', subj_id, run_id)), 'regs');
            load(fullfile(get_mat_dir(true), sprintf('get_visuals_subj%d_run%d.mat', subj_id, run_id)), 'visuals');
            load(fullfile(get_mat_dir(true), sprintf('get_onoff_subj%d_run%d.mat', subj_id, run_id)), 'onoff');

            % extract event onsets
            for i = 1:length(regs_fields)
                field = regs_fields{i};
                onsets(run_id).(field) = regs.timestamps(logical(regs.(field)));
            end
            for i = 1:length(visuals_fields)
                field = visuals_fields{i};
                onsets(run_id).(field) = regs.timestamps(logical(visuals.(field))); % use regs.timestamps, not visuals.timestamps, to cross-reference events (they are slightly off)
            end
            for i = 1:length(onoff_fields)
                field = onoff_fields{i};
                onsets(run_id).(field) = onoff.(field);
            end
        end

        toc

        if strcmp(what, 'GP')
            disp('extracting predicted BOLD time course from EMPA GP results');
            load(fullfile(get_mat_dir(false), sprintf('fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=EMPA_theory_nsamples=100_project=1_norm=1_fast=1_saveYhat=1.mat', subj_id)), 'Y_hat');
        elseif strcmp(what, 'GP_DQN')
            disp('extracting predicted BOLD time course from DQN GP results');
            load(fullfile(get_mat_dir(false), sprintf('fit_gp_CV_HRR_subj=%d_us=1_glm=1_mask=mask_model=DQN_all_nsamples=100_project=1_norm=1_fast=1_saveYhat=1.mat', subj_id)), 'Y_hat');
        end

        disp('extracting BOLD time courses');

        % loop over masks
        for m = 1:length(mask_filenames)
            mask_filename = mask_filenames{m};
            [~, mask_name{m}, ~] = fileparts(mask_filename);
            disp(mask_name{m});

            % initialize commutative timecourses and counts
            % we divide them in the end to get the PETH for a given subject
            % notice that we maintain a separate count for each TR/bin
            % this accounts for "TRs" outside of the run
            for i = 1:length(fields)
                field = fields{i};
                activations(m).(field)(s,:) = zeros(1, length(PETH_dTRs));
                counts(m).(field)(s,:) = zeros(1, length(PETH_dTRs));
            end

            % extract ROI mask
            [mask_format, mask, Vmask] = get_mask_format_helper(mask_filename);

            % get BOLD time course from ROI
            [Y, K, W, R, Y_run_id] = load_BOLD(EXPT, glmodel, subj_id, mask, Vmask);
            % dummy BOLD, for local testing TODO disable
            %Y = rand(EXPT.nTRs * 6, 10);
            %Y_run_id = [ones(EXPT.nTRs, 1) * 1, ...
            %            ones(EXPT.nTRs, 1) * 2, ...
            %            ones(EXPT.nTRs, 1) * 3, ...
            %            ones(EXPT.nTRs, 1) * 4, ...
            %            ones(EXPT.nTRs, 1) * 5, ...
            %            ones(EXPT.nTRs, 1) * 6];
            tic

            % loop over runs
            for SPM_run_id = 1:length(EXPT.subject(subj_id).functional)

                % get whatever we will be plotting on the histograms
                switch what
                    case 'BOLD'
                        % get BOLD time course given run
                        Y_run = nanmean(Y(Y_run_id == SPM_run_id, :), 2);
                        assert(all(size(Y_run) == [EXPT.nTRs, 1]));
                    case {'GP', 'GP_DQN'}
                        % get predicted BOLD and BOLD, do not average across voxels
                        Y_run = Y(Y_run_id == SPM_run_id, :);
                        Y_hat_run = Y_hat(Y_run_id == SPM_run_id, mask(whole_brain_mask));
                        % get correlation across voxels at each time point
                        r_run = nan(EXPT.nTRs, 1);
                        for t = 1:EXPT.nTRs
                            r_run(t) = corr(Y_run(t, :)', Y_hat_run(t, :)');
                        end
                        z_run = atanh(r_run);
                    otherwise
                        assert(false);
                end

                run_id = get_behavioral_run_id(subj_id, SPM_run_id);

                % get BOLD around event onsets
                % loop over event types
                for i = 1:length(fields)
                    field = fields{i};
                    event_onsets = onsets(run_id).(field);

                    % loop over individual events
                    for j = 1:length(event_onsets)
                        % which TRs to plot
                        event_TR = round(event_onsets(j) / EXPT.TR); % event onset
                        TRs = event_TR + PETH_dTRs;
                        valid_TRs = TRs >= 1 & TRs <= EXPT.nTRs;

                        % baseline
                        switch what
                            case 'BOLD'
                                Y_baseline = nanmean(Y_run(event_TR + baseline_dTRs));
                            case {'GP', 'GP_DQN'}
                                z_baseline = nanmean(z_run(event_TR + baseline_dTRs));
                            otherwise
                                assert(false);
                        end

                        % accumulate peri-event timecourses; we average at the end (per subject)
                        switch what
                            case 'BOLD'
                                activations(m).(field)(s,valid_TRs) = activations(m).(field)(s,valid_TRs) + (Y_run(TRs(valid_TRs))' - Y_baseline);
                            case {'GP', 'GP_DQN'}
                                activations(m).(field)(s,valid_TRs) = activations(m).(field)(s,valid_TRs) + (z_run(TRs(valid_TRs))' - z_baseline);
                            otherwise
                                assert(false);
                        end
                        counts(m).(field)(s,:) = counts(m).(field)(s,:) + valid_TRs;
                    end
                end
            end

            toc

            % compute PETH
            for i = 1:length(fields)
                field = fields{i};
                activations(m).(field)(s,:) = activations(m).(field)(s,:) ./ counts(m).(field)(s,:);
            end

        end % loop over masks
    end % loop over subjects

    save(filename);

    disp('done');
