function betas_to_scores(glmodel, contrast, Num, sphere, regressor_name)
    % correlate betas to subject scores

    EXPT = vgdl_expt();

    [~, subjdirs, goodRuns, goodSubjects] = vgdl_getSubjectsDirsAndRuns();
    subj_ids = [1:1:32];

    %sphere = 4; % sphere radius in mm
    %Num = 3; % # peaks per ROI
    %glmodel = 21;
    %contrast = 'theory_change_flag';
    %regressor_name = 'theory_change_flag'

    % spherical mask around top ROI from contrast
    [mask_filenames, regions] = get_masks_from_contrast(glmodel, contrast, true, [], Num, sphere);
    filename = sprintf('mat/betas_to_scores_glm=%d_con=%s_Num=%d_sphere=%.1fmm_reg=%s.mat', glmodel, contrast, Num, sphere, regressor_name);
    disp(filename);

    % which runs 
    score_run_ids = [5 6];
    beta_run_ids = [1 2 3 4 ];

    for m = 1:length(mask_filenames)
        mask_filename = mask_filenames{m};
        [~, mask_name{m}, ~] = fileparts(mask_filename);
        disp(mask_name{m});

        scores{m} = [];
        wins{m} = [];
        success_rates{m} = [];
        betas{m} = [];

        for s = 1:length(subj_ids)
            subj_id = subj_ids(s);
            fprintf('Subject %d\n', subj_id);

            % get subject score from last two runs
            [instance_scores, instance_wins, instant_success_rates, game_names] = get_instance_scores(subj_id, score_run_ids, true);
            score = mean(instance_scores);
            win = mean(instance_wins);
            success_rate = mean(instant_success_rates);

            % get betas from first four runs
            %B = ccnl_get_beta_series(EXPT, glmodel, subj_id, regressor_name, mask_filename);
            B = rand(sum(goodRuns{subj_id}), 124); % for debugging

            assert(size(B, 1) == sum(goodRuns{subj_id}));
            SPM_beta_run_ids = get_SPM_run_ids(subj_id, beta_run_ids);
            SPM_beta_run_ids = SPM_beta_run_ids(~isnan(SPM_beta_run_ids));

            SPM_beta_run_ids
            if isempty(SPM_beta_run_ids)
                fprintf('No SPM_beta_run_ids for subject %d\n', subj_id);
                continue; % no runs...
            end

            B = B(SPM_beta_run_ids, :);
            beta = mean(B(:));
            %beta = max(mean(B, 1));

            % append betas and scores
            scores{m} = [scores{m}, score];
            wins{m} = [wins{m}, win];
            success_rates{m} = [success_rates{m}, success_rate];
            betas{m} = [betas{m}, beta];
        end

		[r,p] = corr(scores{m}', betas{m}', 'type', 'pearson');
		ps(m,1) = p;
		rs(m,1) = r;

		[r,p] = corr(wins{m}', betas{m}', 'type', 'pearson');
		ps(m,2) = p;
		rs(m,2) = r;

		[r,p] = corr(success_rates{m}', betas{m}', 'type', 'pearson');
		ps(m,3) = p;
		rs(m,3) = r;
    end

    p_uncorr = ps;
    p_corr = 1 - (1 - ps) .^ size(ps, 1);
    r = rs;

    save(filename);

    table(mask_filenames', p_uncorr, p_corr, r)
