function [game_names_ordered, sprite_types_ordered] = get_game_names_ordered(subj_id)
    if subj_id <= 11
        game_names_ordered = {'vgfmri3_chase','vgfmri3_helper','vgfmri3_bait','vgfmri3_lemmings','vgfmri3_plaqueAttack','vgfmri3_zelda'};
        sprite_types_ordered = {
            { 'wall', 'carcass', 'angry', 'scared' }, ... % chase 
			{ 'wall', 'box1', 'box2', 'box3', 'forcefield', 'chaser1', 'chaser2', 'missile1', 'missile2' }, ... # helper
			{ 'wall', 'hole', 'box', 'key', 'goal', 'mushroom' }, ... # bait
			{ 'wall', 'goal', 'entrance', 'hole', 'lemming', 'shovel' }, ... # lemmings
			{ 'wall', 'flour', 'hotdog', 'hotdoghole', 'burger', 'burgerhole', 'fullMolarInf' ,'fullMolarSup', 'deadMolarInf', 'deadMolarSup' }, ... # plaqueAttack
			{ 'wall', 'goal', 'key', 'monsterQuick', 'monsterNormal', 'monsterSlow', 'sword' } ... # zelda
        };
    else
        game_names_ordered = {'vgfmri4_chase', 'vgfmri4_helper', 'vgfmri4_bait', 'vgfmri4_lemmings', 'vgfmri4_avoidgeorge', 'vgfmri4_zelda'};
        sprite_types_ordered = {
            { 'wall', 'carcass', 'angry', 'scared' }, ... % chase 
			{ 'wall', 'box1', 'box2', 'box3', 'forcefield', 'chaser1', 'chaser2', 'missile1', 'missile2' }, ... # helper
			{ 'wall', 'hole', 'box', 'key', 'goal', 'mushroom' }, ... # bait
			{ 'wall', 'goal', 'entrance', 'hole', 'lemming', 'shovel' }, ... # lemmings
			{ 'wall', 'quiet', 'george', 'cigarette' }, ... # avoidgeorge
			{ 'wall', 'goal', 'key', 'monsterQuick', 'monsterNormal', 'monsterSlow', 'sword' } ... # zelda
		}
    end

