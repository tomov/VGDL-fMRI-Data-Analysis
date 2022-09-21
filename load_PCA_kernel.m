function [ker, features] = load_PCA_kernel(subj_id, which_run_ids, normalize)
    filename = sprintf('PCA_subject_kernel_subj=%d_sigma_w=1.000_norm=%d.mat', subj_id, normalize);
    filename = fullfile(get_mat_dir(), filename);
    filename
    load(filename, 'projection_kernel', 'r_id', 'projection_Xx');

    ker = projection_kernel;
    features = projection_Xx;
    assert(length(r_id) == size(ker, 1));
    assert(length(r_id) == size(features, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
    features = features(which_TRs, :);
end

