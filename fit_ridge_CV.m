function fit_ridge_CV(subj, use_smooth, glmodel, mask, model_name, what, subsample_only, project, save_Y_hat)
    % copied from fit_gp_CV.m and decode_gp_CV.m

    %{
    clear all;
    subj = 1;
    use_smooth = true;
    glmodel = 1;
    %mask = 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii';
    %mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
    %mask = 'masks/ROI_x=16_y=-94_z=22_1voxels_Sphere1.nii';
    mask = 'masks/N_Acc.nii';
    %mask = 'masks/mask_batchsize=1000_batch=2.nii';
    %mask = 'masks/mask.nii';
    model_name = 'EMPA'
    what = 'theory';
    fast = true;
    project = true;
    debug = false;
    save_Y_hat = true;
    subsample_only = false;
    %}

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

    %what = 'theory';

    %assert(isequal(what, 'theory'));

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end

    if ~exist('save_Y_hat', 'var')
        save_Y_hat = false;
    end

    [~,maskname,~] = fileparts(mask);
    filename = sprintf('fit_ridge_CV_HRR_subj=%d_us=%d_glm=%d_mask=%s_model=%s_%s_subsample=%d_project=%d_saveYhat=%d.mat', subj, use_smooth, glmodel, maskname, model_name, what, subsample_only, project, save_Y_hat);
    filename = fullfile(get_mat_dir(), filename);
    filename


    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);


    % load BOLD
    %
    fprintf('loading BOLD for subj %d\n', subj);
    tic
    [Y, K, W, R, SPM_run_id] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
    run_id = get_behavioral_run_id(subj, SPM_run_id)';
    toc


    % create kernel and HRR regressors from theory id sequence
    %
    fprintf('loading HRRs/embeddings for subj %d\n', subj);
    tic


    switch model_name
        case 'EMPA'
            [Xx, ker] = load_HRR_Xx(subj, unique(run_id), what, subsample_only);

        case 'DQN'
            [Xx] = load_DQN_Xx(subj, unique(run_id), what);

        case 'game'
            [ker, Xx] = load_game_kernel(EXPT, subj); % GLM 1 game id features

        case 'nuisance'
            [ker, Xx] = load_nuisance_kernel(EXPT, subj); % GLM 9 game id features

        case 'nuisance'
            [ker, features] = load_nuisance_kernel(EXPT, subj);

        case 'state'
            [ker, features] = load_state_kernel(EXPT, subj); 

        case 'irrelevant'
            [ker, features] = load_irrelevant_kernel(EXPT, subj); 

        otherwise
            assert(false, 'invalid model name')
    end
    toc


    % whiten, filter & project out nuisance regressors
    %
    if project
        Y = R*K*W*Y;
        Xx = R*K*W*Xx;
    end

    % offset by 3 TRs if not convolving with HRF
    if subsample_only
        Y = Y(4:end,:);
        Xx = Xx(1:end-3,:);
    end

    % every couple of runs form a partition
    partition_id = partition_id_from_run_id(run_id);
    if subsample_only
        partition_id = partition_id(4:end,:);
    end
    assert(size(partition_id, 1) == size(Y, 1));
    n_partitions = max(partition_id);

    disp('run_id');
    disp(run_id');
    disp('partition_id');
    disp(partition_id');


    lambdas = logspace(-5,5,20);

        
    % precompute pseudoinverses for all folds for all lambdas
    % 
    fprintf('precomputing pseudoinverses\n');
    tic
    for j = 1:length(lambdas) % grid search lambda
        l = lambdas(j);

        for k = 1:n_partitions % validate fold
            train = partition_id ~= k;

            % see test_gp.m and spm_spm.m
            X = Xx(train,:);
            pX_CV{j,k} = pinv(X' * X + l * eye(size(X,2))) * X';
        end

        X = Xx;
        pX{j} = pinv(X' * X + l * eye(size(X,2))) * X';
    end

    if save_Y_hat
        Y_hat_CV = nan(size(Y));
    end

    % fit ridge regression
    %

    fprintf('solving ridge for subj %d, %d voxels\n', subj, size(Y,2));

    lambda = nan(1, size(Y,2)); % lambdas
    R2_CV = nan(n_partitions, size(Y,2)); % R^2
    adjR2_CV = nan(n_partitions, size(Y,2)); % adjusted R^2
    r_CV = nan(n_partitions, size(Y,2)); % Pearson correlation
    MSE_CV = nan(n_partitions, size(Y,2)); % MSE
    SMSE_CV = nan(n_partitions, size(Y,2)); % SMSE

    batch_size = 10000;

    % do it in batches
    %
    for i = 1:batch_size:size(Y,2)
        s = i;
        e = i + batch_size - 1;
        if e > size(Y,2)
            e = size(Y,2)
        end

        fprintf('batch %d-%d\n', s, e);
        tic

        y = Y(:,s:e);

        % k-fold CV to pick lambda only
        % for leave-one-subject-out CV for lambda
        %
        y_hat = nan(size(y));
        mses = nan(length(lambdas), size(y,2));

        % pick lambda using CV
        %
        for j = 1:length(lambdas) % grid search lambda
            l = lambdas(j);

            y_pred = nan(size(y));

            for k = 1:n_partitions % validate fold: for fitting lambda
                validate = partition_id == k;
                train = partition_id ~= k;

                % see test_gp.m and spm_spm.m
                %X = Xx(train,:);
                %beta = pinv(X' * X + l * eye(size(X,2))) * X' * y(train,:);
                beta = pX_CV{j,k} * y(train,:);
                y_hat(validate,:) = Xx(validate,:) * beta;
            end

            mses(j,:) = mean((y - y_hat).^2, 1);
        end

        [~,ix] = min(mses, [], 1);
        lambda(s:e) = lambdas(ix);

        % CV evaluate/predict TODO double dipping lambda
        %
        for j = 1:length(ix) % for each voxel
            y_hat = nan(size(y,1),1);

            for k = 1:n_partitions 
                test = partition_id == k;
                train = partition_id ~= k;

                beta = pX_CV{ix(j),k} * y(train,j);
                y_hat(test,:) = Xx(test,:) * beta;

                [R2_CV(k,s+j-1), adjR2_CV(k,s+j-1)] = calc_R2(y(test,j), y_hat(test,:), 1);

                % Pearson
                r_CV(k,s+j-1) = corr(y_hat(test,:), y(test,j));

                % MSE and SMSE, Sec. 2.5 in Rasmussen
                MSE_CV(k,s+j-1) = immse(y(test,j), y_hat(test,:));
                SMSE_CV(k,s+j-1) = MSE_CV(k,s+j-1) / var(y(test,j));
            end

            if save_Y_hat
                Y_hat_CV(:, j + s - 1) = y_hat;
            end
        end

        %{
        % totally overfit
        %
        for j = 1:length(ix)
            beta = pX{ix(j)} * y(:,j);
            y_hat = Xx * beta;

            [R2(s+j-1), adjR2(s+j-1)] = calc_R2(y(:,j), y_hat, 1);

            % Pearson
            r(s+j-1) = corr(y_hat, y(:,j));

            % MSE and SMSE, Sec. 2.5 in Rasmussen
            MSE(s+j-1) = immse(y(:,j), y_hat);
            SMSE(s+j-1) = MSE(s+j-1) / var(y(:,j));
        end
        %}

        toc
    end

    %{
    y = y(:,end);

    figure;
    hold on;
    plot(y);
    plot(y_hat);
    %}

    filename

    if save_Y_hat
        save(filename, 'lambda', 'R2_CV', 'adjR2_CV', 'r_CV', 'MSE_CV', 'SMSE_CV', ...
                       'subj', 'use_smooth', 'glmodel', 'mask', 'what', 'subsample_only', 'project', ...
                       'lambdas', 'X', 'Xx', 'partition_id', 'run_id', 'SPM_run_id', 'Y', 'Y_hat_CV', ...
        '-v7.3');
    else
        save(filename, 'lambda', 'R2_CV', 'adjR2_CV', 'r_CV', 'MSE_CV', 'SMSE_CV', ...
                       'subj', 'use_smooth', 'glmodel', 'mask', 'what', 'subsample_only', 'project', ...
                       'lambdas', 'X', 'Xx', 'partition_id', 'run_id', 'SPM_run_id', ...
        '-v7.3');
    end

    disp('Done');

end

