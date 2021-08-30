function searchlight_rsa2(subj, use_smooth, glmodel, mask, model_name, what, project, sphere)
    % manually do Searchlight RSA
    % TODO deduplicate with fit_gp_CV.m, neurosynth_rsa2.m

    %subj = 1;
    %use_smooth = true;
    %glmodel = 1;
    %%mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
    %%mask = 'masks/mask_batchsize=1000_batch=2.nii';
    %mask = 'masks/mask.nii';
    %model_name = 'EMPA';
    %what = 'theory';
    %project = false;
    %neural_distance = 'correlation';
    %sphere = 10; % mm

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end

    if ~exist('debug', 'var')
        debug = false;
    end
    if ~exist('project', 'var')
        project = false;
    end
    if ~exist('neural_distance', 'var')
        neural_distance = 'correlation';
    end



    [~,maskname,~] = fileparts(mask);
    filename = sprintf('mat/searchlight_rsa2_subj=%d_us=%d_glm=%d_mask=%s_model=%s_%s_nsamples=100_project=%d_dist=%s_r=%.2fmm.mat', subj, use_smooth, glmodel, maskname, model_name, what, project, neural_distance, sphere);
    filename

    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);

    fprintf('loading BOLD for subj %d\n', subj);
    tic
    [Y, K, W, R, SPM_run_id] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
    run_id = get_behavioral_run_id(subj, SPM_run_id)';
    toc

    fprintf('loading kernel for subj %d\n', subj);
    tic
    switch model_name
        case 'EMPA'
            assert(ismember(what, {'theory', 'sprite', 'interaction', 'termination'}));
            ker = load_HRR_kernel(subj, unique(run_id), what);
        case 'DQN'
            assert(ismember(what, {'conv1', 'conv2', 'conv3', 'linear1', 'linear2'}));
            ker = load_DQN_kernel(subj, unique(run_id), what)
        case 'game'
            ker = load_game_kernel(EXPT, subj); % GLM 1 game id features
        case 'nuisance'
            ker = load_nuisance_kernel(EXPT, subj); % GLM 9 game id features
        otherwise
            assert(false, 'invalid model name')
    end
    toc

    fprintf('Memory usage: %.3f MB\n', monitor_memory_whos);

    % whiten, filter & project out nuisance regressors
    if project
        Y = R*K*W*Y;
        ker = R*K*W*ker*W'*K'*R';
    end

    % every couple of runs form a partition
    partition_id = partition_id_from_run_id(run_id);
    assert(size(partition_id, 1) == size(Y, 1));
    n_partitions = max(partition_id);

    disp('run_id');
    disp(run_id');
    disp('partition_id');
    disp(partition_id');





    % new stuff begins here

    % extract behavioral RDM
    assert(ismember(model_name, {'EMPA', 'DQN'}));
    behavioral_RDM = 1 - ker; % assuming the kernel was generated with Z scored representations

    upper = logical(triu(ones(size(behavioral_RDM)), 1));

    r = sphere / 1.5; % mm -> voxels

    % get Searchlight centers
    idx = find(mask);
    [x, y, z] = ind2sub(size(mask), idx);
    cor = [x y z];

    tic
    for i = 1:size(cor, 1)
        if mod(i, 100) == 0
            i
            toc
            tic
        end
        spherical_mask_filename = sprintf('mat/sphere_temp_r=%.2fmm.mat', sphere);
        spherical_mask = ccnl_create_spherical_mask(cor(i,1), cor(i,2), cor(i,3), r, []);
        assert(all(mask(spherical_mask)));

        Y_roi = Y(:, spherical_mask(mask));
        neural_RDM = squareform(pdist(Y_roi, neural_distance));

        [rho(i), p(i)] = corr(neural_RDM(upper), behavioral_RDM(upper), 'type', 'Spearman');
    end

    filename
    save(filename);
