function betas_to_scores_rois(glmodel, regressor_name, atlas)
    % correlate betas to subject scores
    % TODO dedupe w betas_to_scores.m
    conn = mongo('holy7c22108.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54', 'UserName', 'reader', 'Password', 'parolatamadafaqa')

    EXPT = vgdl_expt();

    [~, subjdirs, goodRuns, goodSubjects] = vgdl_getSubjectsDirsAndRuns();
    subj_ids = [1:1:32];

    %sphere = 4; % sphere radius in mm
    %Num = 3; % # peaks per ROI
    %glmodel = 21;
    %contrast = 'theory_change_flag';
    %regressor_name = 'theory_change_flag'

    filename = fullfile(get_mat_dir(false), sprintf('betas_to_scores_rois_glm=%d_reg=%s_atlas=%s.mat', glmodel, regressor_name, atlas));
    disp(filename);

    [mask_filenames, regions, ~] = get_anatomical_masks(atlas);
    regions = regions';

    % which runs 
    score_run_ids = [1 2 3 4 5 6];
    beta_run_ids = [1 2 3 4 5 6];

    for m = 1:length(mask_filenames)
        mask_filename = mask_filenames{m};
        [~, mask_name{m}, ~] = fileparts(mask_filename);
        disp(mask_name{m});

        scores{m} = [];
        wins{m} = [];
        success_rates{m} = [];
        betas{m} = [];
        ts{m} = [];

        for s = 1:length(subj_ids)
            subj_id = subj_ids(s);
            fprintf('Subject %d\n', subj_id);

            % get subject score from last two runs
            [instance_scores, instance_wins, instant_success_rates, game_names] = get_instance_scores(conn, subj_id, score_run_ids, true);
            score = mean(instance_scores);
            win = mean(instance_wins);
            success_rate = mean(instant_success_rates);

            % get betas from first four runs
            B = ccnl_get_beta_series(EXPT, glmodel, subj_id, regressor_name, mask_filename);
            %B = rand(sum(goodRuns{subj_id}), 124); % for debugging
            assert(size(B, 1) == sum(goodRuns{subj_id}));

            SPM_beta_run_ids = get_SPM_run_ids(subj_id, beta_run_ids);
            SPM_beta_run_ids = SPM_beta_run_ids(~isnan(SPM_beta_run_ids));

            SPM_beta_run_ids
            if isempty(SPM_beta_run_ids)
                fprintf('No SPM_beta_run_ids for subject %d\n', subj_id);
                continue; % no runs...
            end

            B = B(SPM_beta_run_ids, :);

            % get betas for all runs
            %{
            B = ccnl_get_beta(EXPT, glmodel, regressor_name, mask_filename, subj_id);
            assert(size(B, 1) == 1);
            %}

            beta = mean(B(:));
            %beta = max(mean(B, 1));

            % get T statistics
            T = ccnl_get_tmap(EXPT, glmodel, regressor_name, mask_filename, subj_id);
            assert(size(T, 1) == 1);
            tstat = mean(T(:));

            % append betas and scores
            scores{m} = [scores{m}, score];
            wins{m} = [wins{m}, win];
            success_rates{m} = [success_rates{m}, success_rate];
            betas{m} = [betas{m}, beta];
            ts{m} = [ts{m}, tstat];
        end

		[r,p] = corr(scores{m}', ts{m}', 'type', 'spearman');
		ps(m,1) = p;
		rs(m,1) = r;

		[r,p] = corr(wins{m}', ts{m}', 'type', 'spearman');
		ps(m,2) = p;
		rs(m,2) = r;

		[r,p] = corr(success_rates{m}', ts{m}', 'type', 'spearman');
		ps(m,3) = p;
		rs(m,3) = r;
    end

    p_uncorr = ps;
    p_corr = 1 - (1 - ps) .^ size(ps, 1);
    r = rs;

    save(filename);

    table(mask_filenames', p_uncorr, p_corr, r)
