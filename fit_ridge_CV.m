%function fit_ridge_CV(subj, use_smooth, glmodel, mask)

%{
clear all;

% copied from fit_gp_CV.m and decode_gp_CV.m

subj = 1;
use_smooth = true;
glmodel = 21;
%mask = 'masks/ROI_x=48_y=12_z=32_62voxels_Sphere6.nii';
%mask = 'masks/ROI_x=48_y=12_z=32_1voxels_Sphere1.nii';
mask = 'masks/ROI_x=42_y=28_z=26_1voxels_Sphere1.nii';
what = 'theory';

assert(isequal(what, 'theory'));

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end


[~,maskname,~] = fileparts(mask);
filename = sprintf('mat/fit_ridge_CV_HRR_subj=%d_us=%d_glm=%d_mask=%s_%s_nsamples=1.mat', subj, use_smooth, glmodel, maskname, what);
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


        


y = Y(:,1);

%}


lambdas = logspace(-3,3,30);

y_hat = nan(size(y));


%{
% nested CV for unbiased predictions
% momchil: not good, seems to underfit => use lambda from other subjects
%
for k = 1:n_partitions % test fold: evaluate fit, for model comparison e.g. w/ GP
    test = partition_id == k;

    % pick lambda using CV
    %
    for i = 1:length(lambdas) % grid search lambda
        lambda = lambdas(i);

        y_pred = nan(size(y));

        for j = 1:n_partitions % validate fold: for fitting lambda
            if k == j
                continue;
            end

            validate = partition_id == j;
            train = partition_id ~= j & partition_id ~= k;

            % see test_gp.m TODO precompute

            X = Xx(train,:);
            beta = pinv(X' * X + lambda * eye(size(X,2))) * X' * y(train,:);
            
            y_pred(validate,:) = Xx(validate,:) * beta;

        end

        mses(i) = immse(y(~test,:), y_pred(~test,:));
    end

    [~,ix] = min(mses);
    lambda = lambdas(ix);
    fold_lambdas(k) = lambda;

    X = Xx(~test,:);
    beta = pinv(X' * X + lambda * eye(size(X,2))) * X' * y(~test,:);
    y_hat(test,:) = Xx(test,:) * beta;
end

%{
% totally overfit
X = Xx;
beta = pinv(X' * X + lambda * var(diag(X))) * X' * y; % note that ridge regression is not scale invariant, that is, X(:,i) * 2 doesn't make b(i) * 2
y_hat = X * beta;
%}


figure;
hold on;
plot(y);
plot(y_hat);
%}



% k-fold CV to pick lambda only
% for leave-one-subject-out CV for lambda
%

y_hat = nan(size(y));

% pick lambda using CV
%
for i = 1:length(lambdas) % grid search lambda
    lambda = lambdas(i);

    y_pred = nan(size(y));

    for j = 1:n_partitions % validate fold: for fitting lambda
        validate = partition_id == j;
        train = partition_id ~= j;

        % see test_gp.m TODO precompute

        X = Xx(train,:);
        beta = pinv(X' * X + lambda * eye(size(X,2))) * X' * y(train,:);
        
        y_pred(validate,:) = Xx(validate,:) * beta;

    end

    mses(i) = immse(y, y_pred);
end

[~,ix] = min(mses);
lambda = lambdas(ix);
fold_lambdas(k) = lambda;


% totally overfit
X = Xx;
beta = pinv(X' * X + lambda * eye(size(X,2))) * X' * y;
y_hat = X * beta;


figure;
hold on;
plot(y);
plot(y_hat);
