% load DQN regressors
%
function [Xx] = load_DQN_Xx(subj_id, which_run_ids, what)
    assert(ismember(what, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2'}));
    load(sprintf('mat/DQN_subject_Xx_subj=%d_sigma_w=1.000_norm=1.mat', subj_id), 'layer_conv1_output_Xx', 'layer_conv2_output_Xx', 'layer_conv3_output_Xx', 'layer_linear1_output_Xx', 'layer_linear2_output_Xx', 'r_id');

    Xx = eval(['layer_', what, '_output_Xx']);
    
    % subset TRs based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    Xx = Xx(which_TRs, :);
end
