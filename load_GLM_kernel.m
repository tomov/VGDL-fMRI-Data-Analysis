function [ker, features] = load_GLM_kernel(EXPT, glmodel, subj_id, regressor_names, do_norm)
    % create GP kernel from GLM & regressor names
    % note that we merge the different runs, b/c we want CV but SPM does random effects across runs

    modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj_id)]);
    load(fullfile(modeldir,'SPM.mat'));

    features = [];
    for i = 1:length(regressor_names)
        which = contains(SPM.xX.name, regressor_names{i});
        feature = sum(SPM.xX.X(:,which), 2); % merge game boxcars from different runs
        features(:,i) = feature;
    end

    % optionally normalize
    if do_norm
        features = zscore(features, 0, 1);
    end

    % see gen_subject_kernels() in HRR.py
    sigma_w = 1; % TODO fit; matching with sigma_w in HRR.py
    Sigma_w = eye(size(features,2)) * sigma_w; % Sigma_p in Rasmussen, Eq. 2.4
    ker = features * Sigma_w * features'; % K in Rasmussen, Eq. 2.12
end
