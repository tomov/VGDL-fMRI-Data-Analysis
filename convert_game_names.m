function [proper_game_names] = convert_game_names(game_names)

    proper_game_names = game_names;
    proper_game_names(strcmp(game_names, 'vgfmri3_chase')) = {'Chase'};
    proper_game_names(strcmp(game_names, 'vgfmri3_helper')) = {'Helper'};
    proper_game_names(strcmp(game_names, 'vgfmri3_bait')) = {'Bait'};
    proper_game_names(strcmp(game_names, 'vgfmri3_lemmings')) = {'Lemmings'};
    proper_game_names(strcmp(game_names, 'vgfmri3_plaqueAttack')) = {'Plaque Attack'};
    proper_game_names(strcmp(game_names, 'vgfmri3_avoidgeorge')) = {'Avoid George'};
    proper_game_names(strcmp(game_names, 'vgfmri3_zelda')) = {'Zelda'};
