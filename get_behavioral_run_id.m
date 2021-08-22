function run_id = get_behavioral_run_id(subj_id, SPM_run_id)
    % get the behavioral run_id from an SPM run_id
    % there can be discrepancies whenever we omit runs from the BOLD analysis due to motion, scanner crashes, etc.
    % the run_id returned here is used to index runs in the behavioral/modeling data
    % the SPM_run_id is used to index runs in the fMRI data

    [allSubjects, subj_dirs, goodRuns, goodSubjs] = vgdl_getSubjectsDirsAndRuns();

    % skip bad run_ids
    assert(find(subj_id == allSubjects) == subj_id); % not sure if this works if we skip subjects
    run_ids = find(goodRuns{find(subj_id == allSubjects)});
    run_id = run_ids(SPM_run_id);
