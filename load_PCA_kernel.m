function [ker] = load_PCA_kernel(subj_id, which_run_ids, normalize)
    filename = sprintf('PCA_subject_kernel_subj=%d_sigma_w=1.000_norm=%d.mat', subj_id, normalize);
    filename = fullfile(get_mat_dir(), filename);
    filename
    load(filename, 'projection_kernel', 'r_id');

    ker = projection_kernel;
    assert(length(r_id) == size(ker, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
end

