% load DQN regressors
%
function [Xx] = load_DQN_Xx(subj_id, which_run_ids, what)
    assert(ismember(what, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2', 'all'}));

    filename = sprintf('DQN_subject_kernel_subj=%d_sigma_w=1.000_norm=1.mat', subj_id);
    filename = fullfile(get_mat_dir(), filename);
    load(filename, 'layer_conv1_output_Xx', 'layer_conv2_output_Xx', 'layer_conv3_output_Xx', 'layer_linear1_output_Xx', 'layer_linear2_output_Xx', 'r_id');

    if strcmp(what, 'all')
        Xx = [layer_conv1_output_Xx layer_conv2_output_Xx layer_conv3_output_Xx layer_linear1_output_Xx layer_linear2_output_Xx];
    else
        Xx = eval(['layer_', what, '_output_Xx']);
    end
    assert(length(r_id) == size(Xx, 1));
    
    % subset TRs based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    Xx = Xx(which_TRs, :);
end
