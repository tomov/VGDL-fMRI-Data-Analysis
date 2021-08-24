function SPM_run_id = get_SPM_run_ids(subj_id, run_id)
    % get the SPM run id from a behavioral run id
    % inverse of get_behavioral_run_id

    [allSubjects, subj_dirs, goodRuns, goodSubjs] = vgdl_getSubjectsDirsAndRuns();

    SPM_run_ids = cumsum(goodRuns{subj_id});
    bad = ~goodRuns{subj_id}(run_id);

    SPM_run_id(bad) = nan;
    SPM_run_id(~bad) = SPM_run_ids(run_id(~bad));

