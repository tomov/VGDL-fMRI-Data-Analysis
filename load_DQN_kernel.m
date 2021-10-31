function [ker] = load_DQN_kernel(subj_id, which_run_ids, what)
    assert(ismember(what, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2', 'all'}));

    filename = sprintf('DQN_subject_kernel_subj=%d_sigma_w=1.000_norm=1.mat', subj_id);
    filename = fullfile(get_mat_dir(), filename);
    load(filename, 'layer_conv1_output_kernel', 'layer_conv2_output_kernel', 'layer_conv3_output_kernel', 'layer_linear1_output_kernel', 'layer_linear2_output_kernel', 'r_id');

    if strcmp(what, 'all')
        % special case for the concatenation of all laye
        % NOTE: this assumes that there is only a single sample,  and that sigma_w = 1
        Xx = [layer_conv1_output_Xx layer_conv2_output_Xx layer_conv3_output_Xx layer_linear1_output_Xx layer_linear2_output_Xx];
        ker = Xx * Xx';
    else
        ker = eval(['layer_', what, '_output_kernel']);
    end
    assert(length(r_id) == size(ker, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
end

