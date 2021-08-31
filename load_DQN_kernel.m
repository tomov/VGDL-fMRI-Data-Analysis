function [ker] = load_DQN_kernel(subj_id, which_run_ids, what)
    filename = sprintf('mat/DQN_subject_kernel_subj=%d_sigma_w=1.000_norm=1.mat', subj_id);
    load(filename, 'layer_conv1_output_kernel', 'layer_conv2_output_kernel', 'layer_conv3_output_kernel', 'layer_linear1_output_kernel', 'layer_linear2_output_kernel', 'r_id');

    ker = eval(['layer_', what, '_output_kernel']);
    assert(length(r_id) == size(ker, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
end