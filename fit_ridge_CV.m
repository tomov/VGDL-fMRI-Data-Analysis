function fit_ridge_CV(subj, use_smooth, glmodel, mask)
    % copied from fit_gp_CV.m and decode_gp_CV.m


    %{
    clear all;

    subj = 1;
    use_smooth = true;
    glmodel = 9;
    mask = 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii';
    %mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
    %mask = 'masks/ROI_x=42_y=28_z=26_1voxels_Sphere1.nii';
    %}

    what = 'theory';

    assert(isequal(what, 'theory'));

    if use_smooth
        EXPT = vgdl_expt();
    else
        EXPT = vgdl_expt_nosmooth();
    end


    [~,maskname,~] = fileparts(mask);
    filename = sprintf('mat/fit_ridge_CV_HRR_subj=%d_us=%d_glm=%d_mask=%s_%s.mat', subj, use_smooth, glmodel, maskname, what);
    filename


    % load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);


    % create kernel and HRR regressors from theory id sequence
    %

    fprintf('loading HRRs for subj %d\n', subj);
    tic

    load(sprintf('mat/unique_HRR_subject_subj=%d_K=10_N=10_E=0.050_nsamples=1_norm=1.mat', subj), 'theory_HRRs', 'run_id', 'ts', 'theory_id_seq', 'play_key', 'gameStrings', 'unique_theories_filename');
    unique_theory_HRRs = theory_HRRs;
    run_id_frames = run_id';
    ts = ts';

    load('mat/SPM73.mat');

    [theory_kernel, ~, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id_frames, SPM);

    toc


    % load BOLD
    %
    fprintf('loading BOLD for subj %d\n', subj);
    tic
    [Y, K, W, R, run_id] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);
    toc

    % whiten, filter & project out nuisance regressors
    %
    Y = R*K*W*Y;
    Xx = R*K*W*Xx;

    X = Xx;


    % get partitions from RSA 3
    rsa = vgdl_create_rsa(3, subj);
    partition_id = rsa.model(1).partitions;
    assert(size(partition_id, 1) == size(Y, 1));
    n_partitions = max(partition_id);

    %}

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


    % fit ridge regression
    %

    fprintf('solving ridge for subj %d, %d voxels\n', subj, size(Y,2));

    lambda = nan(1, size(Y,2)); % lambdas
    R2_CV = nan(3, size(Y,2)); % R^2
    adjR2_CV = nan(3, size(Y,2)); % adjusted R^2
    r_CV = nan(3, size(Y,2)); % Pearson correlation
    MSE_CV = nan(3, size(Y,2)); % MSE
    SMSE_CV = nan(3, size(Y,2)); % SMSE

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

    y = y(:,end);

    figure;
    hold on;
    plot(y);
    plot(y_hat);

    filename

    save(filename, 'lambda', 'R2_CV', 'adjR2_CV', 'r_CV', 'MSE_CV', 'SMSE_CV', ...
                   'subj', 'use_smooth', 'glmodel', 'mask', 'what', 'lambdas', 'HRRs', 'Xx', ...
    '-v7.3');

    disp('Done');

%end
