

[~,~,~,subjs] = vgdl_getSubjectsDirsAndRuns();

%ccnl_fmri_con(vgdl_expt(), 67, ...
%    {'interaction_or_termination_change_flag'}, ...
%     subjs);
%
%ccnl_fmri_con(vgdl_expt(), 68, ...
%    {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag', 'interaction_change_flag + termination_change_flag', 'sprite_change_flag + interaction_change_flag + termination_change_flag', 'interaction_change_flag - termination_change_flag', 'termination_change_flag - interaction_change_flag', 'interaction_change_flag - sprite_change_flag', 'termination_change_flag - sprite_change_flag'}, ...
%	subjs);
%

%ccnl_fmri_con(vgdl_expt(), 69, ...
%    {'likelihood', 'sum_lik_play', 'surprise'}, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 70, ...
%    {'sum_lik_play', 'n_ts', 'sum_lik_play - n_ts'}, ...
%	subjs);
%

%ccnl_fmri_con(vgdl_expt(), 71, ...
%    {'replan_flag'}, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 72, ...
%    {'frames', 'replan_flag'}, ...
%	subjs);
%

ccnl_fmri_con(vgdl_expt(), 73, ...
    {'frames', 'sum_lik'}, ...
	subjs);

ccnl_fmri_con(vgdl_expt(), 74, ...
    {'frames', 'spriteKL'}, ...
	subjs);

ccnl_fmri_con(vgdl_expt(), 75, ...
    {'interaction_or_termination_change_flag'}, ...
	subjs);



%ccnl_fmri_con(vgdl_expt(), 61, ...
%    {'theory_change_flag'}, ...
%     subjs);
%
%ccnl_fmri_con(vgdl_expt(), 60, ...
%    {'theory_change_flag'}, ...
%     subjs);
%
%ccnl_fmri_con(vgdl_expt(), 59, ...
%    {'theory_change_flag'}, ...
%     subjs);
%
%ccnl_fmri_con(vgdl_expt(), 62, ...
%    {'theory_change_flag'}, ...
%     subjs);
%
%ccnl_fmri_con(vgdl_expt(), 63, ...
%    {'theory_change_flag'}, ...
%     subjs);
%
%ccnl_fmri_con(vgdl_expt(), 64, ...
%    {'theory_change_flag'}, ...
%     subjs);
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

%ccnl_fmri_con(vgdl_expt(), 3, ...
%    {'theory_change_flag'}, ...
%     subjs);

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

%ccnl_fmri_con(vgdl_expt(), 21, ...
%    {'theory_change_flag', ...
%     'vgfmri3_bait', 'vgfmri3_chase', 'vgfmri3_helper', 'vgfmri3_lemmings', 'vgfmri3_plaqueAttack', 'vgfmri3_zelda', ... % games
%     'vgfmri3_chase + vgfmri3_helper + vgfmri3_lemmings + vgfmri3_plaqueAttack - vgfmri3_bait - vgfmri3_zelda', ... % agents - no agents
%     'vgfmri3_lemmings + vgfmri3_plaqueAttack + vgfmri3_zelda - vgfmri3_bait - vgfmri3_chase - vgfmri3_helper', ... % shoot vs no shoot
%     'up', 'down', 'left', 'right', 'spacebar', ... % keyholds
%     'frames', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', ... % frame visuals
%     'block_start', 'block_end', 'instance_start', 'instance_end', 'play_start', 'play_end'}, ... % game start/end
%     subjs);


%ccnl_fmri_con(vgdl_expt(), 26, ...
%    {'frames',     'theory_change_flag'}, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 27, ...
%    {'frames',     'sprite_change_flag'}, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 28, ...
%    {'frames',     'termination_change_flag'}, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 29, ...
%    {'frames',     'newEffects_flag'   }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 30, ...
%    {'frames',     'likelihood'        }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 31, ...
%    {'frames',     'sum_lik_play'           }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 57, ...
%    {'frames',     'surprise'           }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 58, ...
%    {'theory_change_flag',   'surprise_flag', 'theory_change_flag - surprise_flag'           }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 55, ...
%    {'theory_change_flag',   'avatar_collision_flag', 'theory_change_flag - avatar_collision_flag'           }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 56, ...
%    {'theory_change_flag',   'collision_flag', 'theory_change_flag - collision_flag'           }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 32, ...
%    {'frames',     'n_ts'              }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 33, ...
%    {'frames',     'num_effects'       }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 34, ...
%    {'frames',     'R_GG'              }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 35, ...
%    {'frames',     'R_GGs_1', 'R_GGs_2',              }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 36, ...
%    {'frames',     'R_SG'              }, ...
%	subjs);
%




%ccnl_fmri_con(vgdl_expt(), 37, ...
%    {'frames',     'R_SGs_1'              }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 38, ...
%    {'frames',     'interaction_change_flag'}, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 39, ...
%    {'frames',     'S_len'             }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 40, ...
%    {'frames',     'I_len'             }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 41, ...
%    {'frames',     'T_len'             }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 42, ...
%    {'frames',     'Igen_len'          }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 43, ...
%    {'frames',     'Tnov_len'          }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 44, ...
%    {'frames',     'Ip_len'            }, ...
%	subjs);
%

%ccnl_fmri_con(vgdl_expt(), 45, ...
%    {'frames',     'dS_len'             }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 46, ...
%    {'frames',     'dI_len'             }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 47, ...
%    {'frames',     'dT_len'             }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 48, ...
%    {'frames',     'dIgen_len'          }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 49, ...
%    {'frames',     'dTnov_len'          }, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 50, ...
%    {'frames',     'dIp_len'            }, ...
%	subjs);
%
%

%ccnl_fmri_con(vgdl_expt(), 51, ...
%    {'interaction_change_flag'}, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 52, ...
%    {'termination_change_flag'}, ...
%	subjs);
%
%
%ccnl_fmri_con(vgdl_expt(), 53, ...
%    {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag', 'interaction_change_flag + termination_change_flag', 'sprite_change_flag + interaction_change_flag + termination_change_flag', 'interaction_change_flag - termination_change_flag', 'termination_change_flag - interaction_change_flag', 'interaction_change_flag - sprite_change_flag', 'termination_change_flag - sprite_change_flag'}, ...
%	subjs);
%
%ccnl_fmri_con(vgdl_expt(), 54, ...
%    {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag', 'interaction_change_flag + termination_change_flag', 'sprite_change_flag + interaction_change_flag + termination_change_flag', 'interaction_change_flag - termination_change_flag', 'termination_change_flag - interaction_change_flag', 'interaction_change_flag - sprite_change_flag', 'termination_change_flag - sprite_change_flag'}, ...
%	subjs);
%
