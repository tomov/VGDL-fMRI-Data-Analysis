% load kernel with state features
%
function [ker, features] = load_state_kernel(EXPT, subj_id)
    % note sprite^ b/c it's contained in other names
    regressor_names = {'sprites^', 'new_sprites', 'killed_sprites', 'collisions', 'effects', 'sprite_groups', 'effectsByCol', 'win', 'loss', 'ended', 'score', 'dscore'};
    [ker, features] = load_GLM_kernel(EXPT, 107, subj_id, regressor_names, true);
end

