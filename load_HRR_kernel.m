% load HRR kernel
%
function [ker] = load_HRR_kernel(subj_id, which_run_ids, what)
    filename = sprintf('HRR_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1.mat', subj_id);
    % load from the ncf Mount
    filename = fullfile(get_mat_dir(true), filename);
    %filename = sprintf('../../VGDL/py_vgdl/mat/HRR_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=1_sigma_w=1.000_norm=1.mat', subj_id);
    load(filename, 'theory_kernel', 'sprite_kernel', 'interaction_kernel', 'termination_kernel', 'r_id');

    ker = eval([what, '_kernel']);
    assert(length(r_id) == size(ker, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
end


