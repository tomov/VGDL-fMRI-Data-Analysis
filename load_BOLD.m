% load BOLD for subject, project out GLM
%
function [Y, K, W, R, run_id] = load_BOLD(EXPT, glmodel, subj_id, mask, Vmask)
    % load subject data
    % Y = raw BOLD
    % K = filter matrix
    % W = whitening matrix
    % R = residual forming matrix
    %
    % 
    % Y = Xb + Gg + e              (Y = BOLD, X = HRR features; G = nuisance regressors (X in GLM 21), e = noise)
    % KWY = KWXb + KWGg + KWe      (whiten & filter w.r.t. G, as in GLM 21)
    % R = (I - KWG (KWG)^+)        (residual forming matrix w.r.t. G from GLM 21)
    % RKWY = RKWXb + 0 + RKWe      (project out nuisance regressors)
    % but...use GP regression instead
    %
    % TODO RKWe not gaussian, also maybe correlated => need to whiten again after GP fit

    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj_id)]);
    load(fullfile(modeldir,'SPM.mat'));
    assert(isempty(Vmask) || isequal(SPM.xY.VY(1).dim, Vmask.dim), 'Different dimensions between mask and activations');
    assert(ndims(mask) < 3 || isequal(SPM.Vbeta(1).dim, size(mask)), 'Different dimensions between mask and betas');

    % extract data and design matrix from confound GLM
    %

    Y = spm_data_read(SPM.xY.VY, find(mask)); % BOLD data

    X = SPM.xX.X; % original design matrix
    K = SPM.xX.K; % high-pass filter
    W = SPM.xX.W; % whitening matrix
    KWX = SPM.xX.xKXs.X; % high-pass filtered & whitened design matrix

    R = spm_sp('r',SPM.xX.xKXs); % residual forming matrix R = I - X * pinv(X)

    KWY = spm_filter(K,W*Y); % high-pass filtered & whitened data

    % convert filter K to matrix form
    % see spm_filter.m
    for s = 1:length(K)
        I(K(s).row,K(s).row) = eye(length(K(s).row));
        X0X0(K(s).row, K(s).row) = K(s).X0*K(s).X0';
    end
    K = I - X0X0; % high-pass filter matrix

    assert(immse(K*W*X, KWX) < 1e-15);
    assert(immse(K*W*Y, KWY) < 1e-15);
    assert(immse(R, eye(size(X,1)) - SPM.xX.xKXs.u*SPM.xX.xKXs.u') < 1e-15);

    for r = 1:length(SPM.Sess)
        run_id(SPM.Sess(r).row,:) = r;
    end
end

