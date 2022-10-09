% load HRR kernel
%
function [ker, features] = load_HRR_states_kernel(subj_id, which_run_ids, normalize)
    filename = sprintf('HRR_states_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=%d.mat', subj_id, normalize)
    filename = fullfile(get_mat_dir(false), filename);
    load(filename, 'state_kernel', 'r_id', 'state_Xx');

    ker = state_kernel;
    features = state_Xx;
    assert(length(r_id) == size(ker, 1));
    assert(length(r_id) == size(features, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
    features = features(which_TRs, :);
end


