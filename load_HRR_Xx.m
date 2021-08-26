% load HRR regressors
%
function [Xx, ker] = load_HRR_Xx(subj_id, which_run_ids, what, subsample_only)
    assert(ismember(what, {'theory', 'sprite', 'interaction', 'termination'}));

    load('mat/SPM73.mat');

    load(sprintf('mat/unique_HRR_subject_subj=%d_K=30_N=100_E=0.050_nsamples=1_norm=1.mat', subj_id), 'theory_HRRs', 'sprite_HRRs', 'interaction_HRRs', 'termination_HRRs', 'run_id', 'ts', 'theory_id_seq', 'play_key', 'gameStrings', 'unique_theories_filename');
    %load(sprintf('/Volumes/fMRI-2/VGDL_rc_mat/unique_HRR_subject_subj=%d_K=30_N=100_E=0.050_nsamples=1_norm=1.mat', subj_id), 'theory_HRRs', 'sprite_HRRs', 'interaction_HRRs', 'termination_HRRs', 'run_id', 'ts', 'theory_id_seq', 'play_key', 'gameStrings', 'unique_theories_filename');

    unique_HRRs = eval([what, '_HRRs']);
    run_id_frames = run_id';
    ts = ts';
    [ker, ~, HRRs, Xx, r_id] = gen_kernel_from_theory_id_seq(unique_HRRs, theory_id_seq, ts, run_id_frames, SPM, subsample_only);

    % subset TRs based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
    Xx = Xx(which_TRs, :);
end

