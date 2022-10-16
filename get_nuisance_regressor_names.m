function regressor_names = get_nuisance_regressor_names(subj_id)
    regressor_names = get_game_names_ordered(subj_id);
    % note sprite^ b/c it's contained in other names
    regressor_names = [regressor_names, {'right',     'up',     'down',     'spacebar',     'left',     'new_sprites',     'killed_sprites',     'sprites^',     'non_walls',     'avatar_moved',     'moved',     'movable',     'collisions',     'effects',     'sprite_groups',     'changed',     'avatar_collision_flag',     'effectsByCol',     'block_start',     'block_end',     'instance_start',     'instance_end',     'play_start',     'play_end'}];
end
