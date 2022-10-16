% load nuisance regressors from GLM 9 
%
function [ker, features] = load_nuisance_kernel(EXPT, subj_id, normalize)
    regressor_names = get_nuisance_regressor_names(subj_id);
    [ker, features] = load_GLM_kernel(EXPT, 101, subj_id, regressor_names, true, normalize);
end

