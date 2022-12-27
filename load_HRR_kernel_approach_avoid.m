% load HRR kernel
%
function [ker, features] = load_HRR_kernel_approach_avoid(subj_id, which_run_ids, what, normalize, concat, novelty)
    filename = sprintf('HRR_approach_avoid_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=%d_concat=%d_novelty=%d.mat', subj_id, normalize, concat, novelty)
    filename = fullfile(get_mat_dir(2), filename);
    load(filename, 'theory_kernel', 'sprite_kernel', 'interaction_kernel', 'termination_kernel', 'novelty_kernel', 'r_id', 'theory_Xx', 'sprite_Xx', 'interaction_Xx', 'termination_Xx', 'novelty_Xx');

    fprintf('loading %s_kernel from %s\n', what, filename);

    ker = eval([what, '_kernel']);
    features = eval([what, '_Xx']);
    assert(length(r_id) == size(ker, 1));
    assert(length(r_id) == size(features, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
    features = features(which_TRs, :);
end


