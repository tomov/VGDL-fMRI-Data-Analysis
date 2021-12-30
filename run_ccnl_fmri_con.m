function run_ccnl_fmri_con(glmodels, subjs)

if ~exist('subjs', 'var')
    [subjs,~,~,~] = vgdl_getSubjectsDirsAndRuns();
end

for i = 1:length(glmodels)
    glmodel = glmodels(i)

    switch (glmodel)

        case 3
            ccnl_fmri_con(vgdl_expt(), 3, ...
                {'theory_change_flag'}, ...
                 subjs);

        case {21}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'theory_change_flag', ...
                 'up', 'down', 'left', 'right', 'spacebar', ... % keyholds
                 'frames', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', ... % frame visuals
                 'block_start', 'block_end', 'instance_start', 'instance_end', 'play_start', 'play_end'}, ... % game start/end
                 subjs);

        case {102}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'theory_change_flag', ...
                 'up', 'down', 'left', 'right', 'spacebar', ... % keyholds
                 'frames', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', ... % frame visuals
                 'play_start', 'play_end'}, ... % game start/end
                 subjs);

        case {9}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                { ...
                 'up', 'down', 'left', 'right', 'spacebar', ... % keyholds
                 'frames', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', ... % frame visuals
                 'block_start', 'block_end', 'instance_start', 'instance_end', 'play_start', 'play_end'}, ... % game start/end
                 subjs);

        case {101}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                { ...
                 'up', 'down', 'left', 'right', 'spacebar', ... % keyholds
                 'frames', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', ... % frame visuals
                 'play_start', 'play_end'}, ... % game start/end
                 subjs);

        case 67
            ccnl_fmri_con(vgdl_expt(), 67, ...
                {'interaction_or_termination_change_flag'}, ...
                 subjs);
            
        case {68, 106}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag', 'interaction_change_flag + termination_change_flag', 'sprite_change_flag + interaction_change_flag + termination_change_flag', 'interaction_change_flag - termination_change_flag', 'termination_change_flag - interaction_change_flag', 'interaction_change_flag - sprite_change_flag', 'termination_change_flag - sprite_change_flag'}, ...
            	subjs);
            

            %ccnl_fmri_con(vgdl_expt(), 69, ...
            %    {'likelihood', 'sum_lik_play', 'surprise'}, ...
            %	subjs);
            %
            %ccnl_fmri_con(vgdl_expt(), 70, ...
            %    {'sum_lik_play', 'n_ts', 'sum_lik_play - n_ts'}, ...
            %	subjs);
            %

        case 71
            ccnl_fmri_con(vgdl_expt(), 71, ...
                {'replan_flag'}, ...
                subjs);
            
        case 72
            ccnl_fmri_con(vgdl_expt(), 72, ...
                {'frames', 'replan_flag'}, ...
                subjs);
            
        case 81
            ccnl_fmri_con(vgdl_expt(), 81, ...
                {'replan_flag'}, ...
                subjs);
            

            %ccnl_fmri_con(vgdl_expt(), 73, ...
            %    {'frames', 'sum_lik'}, ...
            %	subjs);
            %

        case 74
            ccnl_fmri_con(vgdl_expt(), 74, ...
                {'frames', 'spriteKL'}, ...
                subjs);
            
        case 76
            ccnl_fmri_con(vgdl_expt(), 76, ...
                {'spriteKL'}, ...
                subjs);

        case 90
            ccnl_fmri_con(vgdl_expt(), 90, ...
                {'spriteKL_flag'}, ...
                subjs);
            %ccnl_fmri_con(vgdl_expt(), 75, ...
            %    {'interaction_or_termination_change_flag'}, ...
            %	subjs);
            %
            %ccnl_fmri_con(vgdl_expt(), 84, ...
            %    {'likelihood'}, ...
            %	subjs);



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



        case 26
            ccnl_fmri_con(vgdl_expt(), 26, ...
                {'frames',     'theory_change_flag'}, ...
                subjs);

        case 27
            ccnl_fmri_con(vgdl_expt(), 27, ...
                {'frames',     'sprite_change_flag'}, ...
            	subjs);
            
        case 28
            ccnl_fmri_con(vgdl_expt(), 28, ...
                {'frames',     'termination_change_flag'}, ...
                subjs);
           
        case 29
            ccnl_fmri_con(vgdl_expt(), 29, ...
                {'frames',     'newEffects_flag'   }, ...
            	subjs);
           
        case 30
            ccnl_fmri_con(vgdl_expt(), 30, ...
                {'frames',     'likelihood'        }, ...
            	subjs);
           
        case 31
            ccnl_fmri_con(vgdl_expt(), 31, ...
                {'frames',     'sum_lik_play'           }, ...
            	subjs);
            
        case 79
            ccnl_fmri_con(vgdl_expt(), 79, ...
                {    'surprise_flag'           }, ...
                subjs);

        case 80
            ccnl_fmri_con(vgdl_expt(), 80, ...
                {    'surprise_flag'           }, ...
                subjs);

        case 85
            ccnl_fmri_con(vgdl_expt(), 85, ...
                {    'sprite_change_flag'           }, ...
                subjs);

        case {86, 103, 159}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {    'sprite_change_flag'           }, ...
                subjs);

        case 57
            ccnl_fmri_con(vgdl_expt(), 57, ...
                {'frames',     'surprise'           }, ...
            	subjs);
            
        case 58
            ccnl_fmri_con(vgdl_expt(), 58, ...
                {'theory_change_flag',   'surprise_flag', 'theory_change_flag - surprise_flag'           }, ...
            	subjs);
            
        case 55
            ccnl_fmri_con(vgdl_expt(), 55, ...
                {'theory_change_flag',   'avatar_collision_flag', 'theory_change_flag - avatar_collision_flag'           }, ...
                subjs);

        case 56
            ccnl_fmri_con(vgdl_expt(), 56, ...
                {'theory_change_flag',   'collision_flag', 'theory_change_flag - collision_flag'           }, ...
                subjs);

        case 32
            ccnl_fmri_con(vgdl_expt(), 32, ...
                {'frames',     'n_ts'              }, ...
            	subjs);
            
        case 33
            ccnl_fmri_con(vgdl_expt(), 33, ...
                {'frames',     'num_effects'       }, ...
            	subjs);
            
        case 34
            ccnl_fmri_con(vgdl_expt(), 34, ...
                {'frames',     'R_GG'              }, ...
            	subjs);
            
        case 35
            ccnl_fmri_con(vgdl_expt(), 35, ...
                {'frames',     'R_GGs_1',            }, ...
            	subjs);
            
        case 36
            ccnl_fmri_con(vgdl_expt(), 36, ...
                {'frames',     'R_SG'              }, ...
            	subjs);

        case 37
            ccnl_fmri_con(vgdl_expt(), 37, ...
                {'frames',     'R_SGs_1'              }, ...
            	subjs);
            
        case 38
            ccnl_fmri_con(vgdl_expt(), 38, ...
                {'frames',     'interaction_change_flag'}, ...
                subjs);

        case {51, 104, 160}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'interaction_change_flag'}, ...
                subjs);

        case 82
            ccnl_fmri_con(vgdl_expt(), 82, ...
                {'interaction_change_flag'}, ...
                subjs);
            
        case 87
            ccnl_fmri_con(vgdl_expt(), 87, ...
                {'interaction_change_flag', 'termination_change_flag'}, ...
                subjs);

        case 88
            ccnl_fmri_con(vgdl_expt(), 88, ...
                {'interaction_change_flag', 'termination_change_flag'}, ...
                subjs);

        case 39
            ccnl_fmri_con(vgdl_expt(), 39, ...
                {'frames',     'S_len'             }, ...
            	subjs);
            
        case 40
            ccnl_fmri_con(vgdl_expt(), 40, ...
                {'frames',     'I_len'             }, ...
            	subjs);
            
        case 41
            ccnl_fmri_con(vgdl_expt(), 41, ...
                {'frames',     'T_len'             }, ...
            	subjs);
            
        case 42
            ccnl_fmri_con(vgdl_expt(), 42, ...
                {'frames',     'Igen_len'          }, ...
            	subjs);
            
        case 43
            ccnl_fmri_con(vgdl_expt(), 43, ...
                {'frames',     'Tnov_len'          }, ...
            	subjs);
            
        case 44
            ccnl_fmri_con(vgdl_expt(), 44, ...
                {'frames',     'Ip_len'            }, ...
            	subjs);
            

        case 45
            ccnl_fmri_con(vgdl_expt(), 45, ...
                {'frames',     'dS_len'             }, ...
            	subjs);
            
        case 46
            ccnl_fmri_con(vgdl_expt(), 46, ...
                {'frames',     'dI_len'             }, ...
            	subjs);
            
        case 47
            ccnl_fmri_con(vgdl_expt(), 47, ...
                {'frames',     'dT_len'             }, ...
            	subjs);
            
        case 48
            ccnl_fmri_con(vgdl_expt(), 48, ...
                {'frames',     'dIgen_len'          }, ...
            	subjs);
            
        case 49
            ccnl_fmri_con(vgdl_expt(), 49, ...
                {'frames',     'dTnov_len'          }, ...
            	subjs);
            
        case 50
            ccnl_fmri_con(vgdl_expt(), 50, ...
                {'frames',     'dIp_len'            }, ...
            	subjs);
            
        case {52, 105}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'termination_change_flag'}, ...
                subjs);

        case 83
            ccnl_fmri_con(vgdl_expt(), 83, ...
                {'termination_change_flag'}, ...
                subjs);
            
            
        case 53 
            ccnl_fmri_con(vgdl_expt(), 53, ...
                {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag', 'interaction_change_flag + termination_change_flag', 'sprite_change_flag + interaction_change_flag + termination_change_flag', 'interaction_change_flag - termination_change_flag', 'termination_change_flag - interaction_change_flag', 'interaction_change_flag - sprite_change_flag', 'termination_change_flag - sprite_change_flag'}, ...
                subjs);
           
        case 54
            ccnl_fmri_con(vgdl_expt(), 54, ...
                {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag', 'interaction_change_flag + termination_change_flag', 'sprite_change_flag + interaction_change_flag + termination_change_flag', 'interaction_change_flag - termination_change_flag', 'termination_change_flag - interaction_change_flag', 'interaction_change_flag - sprite_change_flag', 'termination_change_flag - sprite_change_flag'}, ...
                subjs);
            
        case {112, 114}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'exploratory_goals'}, ...
                subjs);

        case {113, 115}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'exploitative_goals'}, ...
                subjs);

        case {123, 124}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'exploitative_goals', 'exploratory_goals', 'exploitative_goals - exploratory_goals', 'exploratory_goals - exploitative_goals'}, ...
                subjs);

        case {116, 117, 118, 119, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'exploit', 'explore', 'exploit - explore', 'explore - exploit'}, ...
                subjs);

        case {149, 150, 152, 153}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'exploit_start', 'exploit_end', 'explore_start', 'explore_end', 'exploit_start - explore_start', 'explore_start - exploit_start', 'exploit_end - explore_end', 'explore_end - exploit_end', 'exploit_start - exploit_end', 'explore_start - explore_end', 'explore_start + exploit_start'}, ...
                subjs);

        case {151}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'exploit_start', 'explore_start', 'exploit_start - explore_start', 'explore_start - exploit_start', 'keypresses', 'subgoal'}, ...
                subjs);

        case {157, 158}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'new_theory_change_flag', ...
                 'up', 'down', 'left', 'right', 'spacebar', ... % keyholds
                 'frames', 'new_sprites', 'killed_sprites', 'sprites', 'non_walls', 'avatar_moved', 'moved', 'movable', 'collisions', 'effects', 'sprite_groups', 'changed', ... % frame visuals
                 'play_start', 'play_end'}, ... % game start/end
                 subjs);

        case {161}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'new_termination_change_flag'}, ...
                subjs);

        case {162}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'falsified_terminations_flag'}, ...
                subjs);

        case {163}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'sprite_change_flag', 'interaction_change_flag', 'new_termination_change_flag'}, ...
                subjs);

        case {164}
            ccnl_fmri_con(vgdl_expt(), glmodel, ...
                {'sprite_change_flag', 'interaction_change_flag', 'falsified_terminations_flag'}, ...
                subjs);

        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end

end
