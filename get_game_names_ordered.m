function game_names_ordered = get_game_names_ordered(subj_id)
    if subj_id <= 11
        game_names_ordered = {'vgfmri3_chase','vgfmri3_helper','vgfmri3_bait','vgfmri3_lemmings','vgfmri3_plaqueAttack','vgfmri3_zelda'};
    else
        game_names_ordered = {'vgfmri4_chase', 'vgfmri4_helper', 'vgfmri4_bait', 'vgfmri4_lemmings', 'vgfmri4_avoidgeorge', 'vgfmri4_zelda'};
    end
