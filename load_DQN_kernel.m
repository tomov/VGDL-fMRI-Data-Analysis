function [ker, features] = load_DQN_kernel(subj_id, which_run_ids, what, normalize, suffix)
    assert(ismember(what, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2', 'all'}));

    %filename = sprintf('DQN_subject_kernel_subj=%d_sigma_w=1.000_norm=%d.mat', subj_id, normalize);
    filename = sprintf('DQN%s_subject_kernel_subj=%d_sigma_w=1.000_norm=%d.mat', suffix, subj_id, normalize); % <--- 25M! in submission
    %filename = sprintf('DQN%s_subject_kernel_subj=%d_sigma_w=1.000_norm=%d_comp=50.mat', suffix, subj_id, normalize);
    filename = fullfile(get_mat_dir(2), 'dqn_layers', filename);

    if strcmp(what, 'all')
        % special case for the concatenation of all layer
        % NOTE: this assumes that there is only a single sample,  and that sigma_w = 1
        load(filename, 'layer_conv1_output_Xx', 'layer_conv2_output_Xx', 'layer_conv3_output_Xx', 'layer_linear1_output_Xx', 'layer_linear2_output_Xx', 'r_id');
        Xx = [layer_conv1_output_Xx layer_conv2_output_Xx layer_conv3_output_Xx layer_linear1_output_Xx layer_linear2_output_Xx];

        ker = Xx * Xx';
        features = Xx;
    else
        load(filename, 'layer_conv1_output_kernel', 'layer_conv2_output_kernel', 'layer_conv3_output_kernel', 'layer_linear1_output_kernel', 'layer_linear2_output_kernel', 'layer_conv1_output_Xx', 'layer_conv2_output_Xx', 'layer_conv3_output_Xx', 'layer_linear1_output_Xx', 'layer_linear2_output_Xx', 'r_id');

        ker = eval(['layer_', what, '_output_kernel']);
        features = eval(['layer_', what, '_output_Xx']);
    end
    assert(length(r_id) == size(ker, 1));
    assert(length(r_id) == size(features, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
    features = features(which_TRs, :);
end

