function fit_gp_CV_noHRF(subj, use_smooth, glmodel, mask, subsample_only)

%{
    clear all;

    subj = 1;
    use_smooth = true;
    glmodel = 9;
    mask = 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii';
    %mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
    %mask = 'masks/ROI_x=42_y=28_z=26_1voxels_Sphere1.nii';
    subsample_only = true;
    %}

    what = 'theory';

    assert(isequal(what, 'theory'));

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end


    [~,maskname,~] = fileparts(mask);
    filename = sprintf('mat/fit_gp_CV_noHRF_HRR_subj=%d_us=%d_glm=%d_mask=%s_subsample=%d_%s.mat', subj, use_smooth, glmodel, maskname, subsample_only, what);
    filename


    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);


    % create kernel and HRR regressors from theory id sequence
    %

    fprintf('loading HRRs for subj %d\n', subj);
    tic

    load(sprintf('mat/unique_HRR_subject_subj=%d_K=10_N=10_E=0.050_nsamples=100_norm=1.mat', subj), 'theory_HRRs', 'run_id', 'ts', 'theory_id_seq', 'play_key', 'gameStrings', 'unique_theories_filename');
    unique_theory_HRRs = theory_HRRs;
    run_id_frames = run_id';
    ts = ts';

    load('mat/SPM73.mat');

    [theory_kernel, ~, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id_frames, SPM, subsample_only);

    % load BOLD
    %
    fprintf('loading BOLD for subj %d\n', subj);
    tic
    [Y, K, W, R, run_id_TRs] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
    toc

    % whiten, filter & project out nuisance regressors
    Y = R*K*W*Y;
    ker = R*K*W*theory_kernel*W'*K'*R';

    % offset by 3 TRs if not convolving with HRF
    if subsample_only
        Y = Y(4:end,:);
        run_id_TRs = run_id_TRs(4:end,:);
        ker = ker(1:end-3,1:end-3);
    end

    % get partitions from RSA 3
    rsa = vgdl_create_rsa(3, subj);
    partition_id = rsa.model(1).partitions;
    if subsample_only
        partition_id = partition_id(4:end,:);
    end
    assert(size(partition_id, 1) == size(Y, 1));
    n_partitions = max(partition_id);

    % init GP stuff
    %
    n = size(ker, 1);
    x = [1:n]';

    % init functions and  hyperparams from kernel
    %meanfun = {@meanDiscrete, n};
    meanfun = @meanConst;
    covfun = {@covDiscrete, n};
    likfun = @likGauss;


    sigmas = 1;

    y = Y;
    [r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj, use_smooth, glmodel, Y, ker, run_id_TRs, x, y, meanfun, covfun, likfun, partition_id);


    save(filename, 'sigmas', 'R2_CV', 'r_CV', 'MSE_CV', 'SMSE_CV', ...
                   'subj', 'use_smooth', 'glmodel', 'mask', 'what', ...
    '-v7.3');

    disp('Done');
