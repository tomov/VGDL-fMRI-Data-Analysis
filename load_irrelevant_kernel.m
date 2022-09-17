% load kernel with irrelevant state features
%
function [ker, features] = load_irrelevant_kernel(EXPT, subj_id, normalize)
    % note sprite^ b/c it's contained in other names
    regressor_names = {'non_walls', 'avatar_moved', 'moved', 'movable', 'changed', 'avatar_collision_flag'};
    [ker, features] = load_GLM_kernel(EXPT, 108, subj_id, regressor_names, true, normalize);
end

