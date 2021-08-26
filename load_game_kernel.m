% load game identity kernel based on GLM 1
%
function [ker, features] = load_game_kernel(EXPT, subj_id)
    game_names = get_game_names_ordered(subj_id);
    [ker, features] = load_GLM_kernel(EXPT, 1, subj_id, game_names, false);
end
