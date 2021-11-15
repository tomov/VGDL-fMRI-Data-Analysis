%function betas_to_scores(glmodel, contrast, Num, sphere, regressor_name)
    % correlate GP statistics to subject scores

    close all;
    clear all;

    EXPT = vgdl_expt();

    [~, subjdirs, goodRuns, goodSubjects] = vgdl_getSubjectsDirsAndRuns();
    subj_ids = [1:1:32];

    out_filename = fullfile(get_mat_dir(false), sprintf('gp_CV_to_scores.mat'));
    out_filename

    agg_filename = fullfile(get_mat_dir(false), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
    agg_filename
    load(agg_filename);

    reg = find(strcmp(regressor_names, 'theory'));
    assert(length(reg) == 1);

    % which runs 
    score_run_ids = [1 2 3 4 5 6];

    for m = 1:nROIs
        roi_names{m}

        scores{m} = [];
        wins{m} = [];
        success_rates{m} = [];
        z{m} = [];
        f{m} = [];

        for s = 1:nsubjects
            subj_id = subjects(s);
            fprintf('Subject %d\n', subj_id);

            % get subject score from last two runs
            [instance_scores, instance_wins, instant_success_rates, game_names] = get_instance_scores([], subj_id, score_run_ids, true);
            score = mean(instance_scores);
            win = mean(instance_wins);
            success_rate = mean(instant_success_rates);

            % append betas and scores
            scores{m} = [scores{m}, score];
            wins{m} = [wins{m}, win];
            success_rates{m} = [success_rates{m}, success_rate];
            z{m} = [z{m}, zs(m,reg,s)];
            f{m} = [f{m}, fs(m,reg,s)];
        end

		[r,p] = corr(scores{m}', f{m}', 'type', 'spearman');
		ps(m,1) = p;
		rr(m,1) = r;

		[r,p] = corr(wins{m}', f{m}', 'type', 'spearman');
		ps(m,2) = p;
		rr(m,2) = r;

		[r,p] = corr(success_rates{m}', f{m}', 'type', 'spearman');
		ps(m,3) = p;
		rr(m,3) = r;

        %{
		[r,p] = corr(scores{m}', z{m}', 'type', 'pearson');
		ps(m,1) = p;
		rr(m,1) = r;

		[r,p] = corr(wins{m}', z{m}', 'type', 'pearson');
		ps(m,2) = p;
		rr(m,2) = r;

		[r,p] = corr(success_rates{m}', z{m}', 'type', 'pearson');
		ps(m,3) = p;
		rr(m,3) = r;
        %}
    end

    p_uncorr = ps;
    p_corr = 1 - (1 - ps) .^ size(ps, 1);
    r = rr;

    out_filename
    save(out_filename);

    table(roi_names', p_uncorr, p_corr, r)
