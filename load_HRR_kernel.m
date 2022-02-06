% load HRR kernel
%
function [ker] = load_HRR_kernel(subj_id, which_run_ids, what, normalize, concat, novelty)
    if concat == 0 && novelty == 1
        % load from the ncf Mount
        %filename = sprintf('HRR_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=%d.mat', subj_id, normalize);
        %filename = fullfile(get_mat_dir(true), filename);
        % load from cannon repro
        filename = sprintf('HRR_cannon_repro_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=%d_concat=0_novelty=1.mat',  subj_id, normalize);
        filename = fullfile(get_mat_dir(0), filename);
        load(filename, 'theory_kernel', 'sprite_kernel', 'interaction_kernel', 'termination_kernel', 'r_id');
    else
        % new stuff
        filename = sprintf('HRR_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=%d_concat=%d_novelty=%d.mat', subj_id, normalize, concat, novelty)
        filename = fullfile(get_mat_dir(false), filename);
        load(filename, 'theory_kernel', 'sprite_kernel', 'interaction_kernel', 'termination_kernel', 'novelty_kernel', 'r_id');
    end

    fprintf('loading %s_kernel from %s\n', what, filename);

    ker = eval([what, '_kernel']);
    assert(length(r_id) == size(ker, 1));

    % subset kernel based on good runs only
    which_TRs = ismember(r_id, which_run_ids);
    ker = ker(which_TRs, which_TRs);
end


