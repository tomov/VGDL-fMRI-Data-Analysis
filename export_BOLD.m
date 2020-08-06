% export raw BOLD for Jiajia and Daphne
%

glmodel = 9;

use_smooth = false;

if use_smooth
    EXPT = vgdl_expt();
    maskfile = 'masks/mask.nii';
    suffix = 'smooth';
else
    EXPT = vgdl_expt_nosmooth();
    maskfile = 'masks/mask_nosmooth.nii';
    suffix = 'nosmooth';
end

subjects = 1:length(EXPT.subject);

[mask, Vmask] = ccnl_load_mask(maskfile);

for s = 1:length(subjects)


    subj = subjects(s);
    fprintf('    subj %d\n', subj);

    % extract game, run, level, block sequences, theory change flags, etc.

    % TODO dedupe with decode_gp_CV

    timestamps = [];
    collision_flags = [];
    avatar_collision_flags = [];
    interaction_change_flags = [];
    termination_change_flags = [];
    theory_change_flags = [];
    game_ids = [];
    block_ids = [];
    run_ids = [];
    level_ids = [];

    for r = 1:6
        query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, r); % in python we index runs from 0 (but not subjects) 

        try
            run = find(conn, 'runs', 'query', query)
            assert(length(run) == 1);

            regs = get_regressors(subj_id, run, conn, true);

            [~, visuals] = get_visuals(subj_id, run, conn, true);

        catch e
            e
            fprintf('loading from cache...\n');
            load(sprintf('mat/get_regressors_subj%d_run%d.mat', subj_id, r));
            load(sprintf('mat/get_visuals_subj%d_run%d.mat', subj_id, r));
        end
            
        assert(length(visuals.durations) == length(regs.timestamps));

        timestamps = [timestamps; regs.timestamps];

        tc = regs.theory_change_flag;
        tc = logical(tc);
        assert(all(regs.timestamps(tc) == regs.theory_change_flag_onsets));

        game_ids = [game_ids; regs.game_ids];
        block_ids = [block_ids; regs.block_ids];
        level_ids = [level_ids; regs.level_ids];
        run_ids = [run_ids; ones(length(regs.theory_change_flag),1) * r];

        theory_change_flags = [theory_change_flags; regs.theory_change_flag];
        interaction_change_flags = [interaction_change_flags; regs.interaction_change_flag];
        termination_change_flags = [termination_change_flags; regs.termination_change_flag];

    end
    theory_change_flags = logical(theory_change_flags);
    interaction_change_flags = logical(interaction_change_flags);
    termination_change_flags = logical(termination_change_flags);


    n = EXPT.nTRs * length(unique(run_ids));
    assert(n == 1698);

    % for each frame, what's the corresponding TR
    TRs = round(timestamps / 2) + (run_ids - 1) * EXPT.nTRs;

    % init TR labels
    block = nan(n, 1);
    run = nan(n, 1);
    level = nan(n, 1);
    game = nan(n, 1);
    theory_change_flag_cnt = zeros(n, 1);
    interaction_change_flag_cnt = zeros(n, 1);
    termination_change_flag_cnt = zeros(n, 1);

    % get continuous labels
    for tr = 1:n
        block(tr) = round(mean(block_ids(TRs == tr)));
        run(tr) = round(mean(run_ids(TRs == tr)));
        level(tr) = round(mean(level_ids(TRs == tr)));
        game(tr) = round(mean(game_ids(TRs == tr)));
    end

    % actual levels
    % TODO hardcoded get 
    actual_level = nan(n,1);
    assert(length(unique(game(~isnan(game)))) == 6);
    for g = 1:6
        rs = unique(run(game == g));
        assert(length(rs) == 3);
        for i = 1:3
            r = rs(i);
            actual_level(game == g & run == r) = level(game == g & run == r) + 1 + (i-1) * 3;
        end
    end
    instance = level;
    level = actual_level;

    % get impulse labels
    frame_idx = find(theory_change_flags);
    for i = 1:length(frame_idx)
        tr = TRs(frame_idx(i));
        theory_change_flag_cnt(tr) = theory_change_flag_cnt(tr) + 1;
    end
    theory_change_flag = theory_change_flag_cnt > 0;

    frame_idx = find(interaction_change_flags);
    for i = 1:length(frame_idx)
        tr = TRs(frame_idx(i));
        interaction_change_flag_cnt(tr) = interaction_change_flag_cnt(tr) + 1;
    end
    interaction_change_flag = interaction_change_flag_cnt > 0;

    frame_idx = find(termination_change_flags);
    for i = 1:length(frame_idx)
        tr = TRs(frame_idx(i));
        termination_change_flag_cnt(tr) = termination_change_flag_cnt(tr) + 1;
    end
    termination_change_flag = termination_change_flag_cnt > 0;

    game_names_ordered = regs.game_names_ordered;
    game_name_to_id = regs.game_name_to_id;

    % BOLD
    [B, runs] = ccnl_get_activations(EXPT, glmodel, mask, subj, false, false);
    filename = sprintf('mat/BOLD_subj%d_%s.mat', glmodel, subj, suffix);

    filename
    save(filename, 'B', 'runs', 'mask', 'Vmask', 'block', 'run', 'instance', 'level', 'game', 'theory_change_flag_cnt', 'theory_change_flag', 'interaction_change_flag_cnt', 'interaction_change_flag', 'termination_change_flag_cnt', 'termination_change_flag', '-v7.3');


    % residuals
    [R] = ccnl_get_residuals(EXPT, glmodel, mask, subj);
    R = R{1};
    filename = sprintf('mat/residuals_glm%d_subj%d_%s.mat', glmodel, subj, suffix);

    filename
    save(filename, 'R', 'runs', 'mask', 'Vmask', 'block', 'run', 'instance', 'level', 'game', 'theory_change_flag_cnt', 'theory_change_flag', 'interaction_change_flag_cnt', 'interaction_change_flag', 'termination_change_flag_cnt', 'termination_change_flag', '-v7.3');


end

