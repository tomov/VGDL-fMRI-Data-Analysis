

[~,~,~,subjs] = vgdl_getSubjectsDirsAndRuns();

%ccnl_fmri_con(vgdl_expt(), 1, ...
%    {'vgfmri3_bait', 'vgfmri3_chase', 'vgfmri3_helper', 'vgfmri3_lemmings', 'vgfmri3_plaqueAttack', 'vgfmri3_zelda', ...
%     'vgfmri3_chase + vgfmri3_helper + vgfmri3_lemmings + vgfmri3_plaqueAttack - vgfmri3_bait - vgfmri3_zelda', ... % agents - no agents
%     'vgfmri3_lemmings + vgfmri3_plaqueAttack + vgfmri3_zelda - vgfmri3_bait - vgfmri3_chase - vgfmri3_helper'}, ... % shoot vs no shoot
%     subjs);

%ccnl_fmri_con(vgdl_expt(), 5, {'up', 'down', 'left', 'right', 'spacebar'}, subjs);
%ccnl_fmri_con(vgdl_expt(), 6, {'keypresses_up', 'keypresses_down', 'keypresses_left', 'keypresses_right', 'keypresses_spacebar'}, subjs);

%ccnl_fmri_con(vgdl_expt(), 7, ...
%    {'frames', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed'}, ...
%     subjs);

ccnl_fmri_con(vgdl_expt(), 3, ...
    {'theory_change_flag'}, ...
     subjs);

%ccnl_fmri_con(vgdl_expt(), 8, ...
%    {'block_start', 'block_end', 'instance_start', 'instance_end', 'play_start', 'play_end'}, ...
%     subjs);

%ccnl_fmri_con(vgdl_expt(), 9, ...
%    {'vgfmri3_bait', 'vgfmri3_chase', 'vgfmri3_helper', 'vgfmri3_lemmings', 'vgfmri3_plaqueAttack', 'vgfmri3_zelda', ... % games
%     'vgfmri3_chase + vgfmri3_helper + vgfmri3_lemmings + vgfmri3_plaqueAttack - vgfmri3_bait - vgfmri3_zelda', ... % agents - no agents
%     'vgfmri3_lemmings + vgfmri3_plaqueAttack + vgfmri3_zelda - vgfmri3_bait - vgfmri3_chase - vgfmri3_helper', ... % shoot vs no shoot
%     'up', 'down', 'left', 'right', 'spacebar', ... % keyholds
%     'frames', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', ... % frame visuals
%     'block_start', 'block_end', 'instance_start', 'instance_end', 'play_start', 'play_end'}, ... % game start/end
%     subjs);
%
