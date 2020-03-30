

[~,~,~,subjs] = vgdl_getSubjectsDirsAndRuns();

%ccnl_fmri_con(vgdl_expt(), 1, ...
%    {'vgfmri3_bait', 'vgfmri3_chase', 'vgfmri3_helper', 'vgfmri3_lemmings', 'vgfmri3_plaqueAttack', 'vgfmri3_zelda', ...
%     'vgfmri3_chase + vgfmri3_helper + vgfmri3_lemmings + vgfmri3_plaqueAttack - vgfmri3_bait - vgfmri3_zelda', ... % agents - no agents
%     'vgfmri3_lemmings + vgfmri3_plaqueAttack + vgfmri3_zelda - vgfmri3_bait - vgfmri3_chase - vgfmri3_helper'}, ... % shoot vs no shoot
%     subjs);

%ccnl_fmri_con(vgdl_expt(), 5, {'up', 'down', 'left', 'right', 'spacebar'}, subjs);
%ccnl_fmri_con(vgdl_expt(), 6, {'keypresses_up', 'keypresses_down', 'keypresses_left', 'keypresses_right', 'keypresses_spacebar'}, subjs);

ccnl_fmri_con(vgdl_expt(), 7, ...
    {'vgfmri3_bait', 'vgfmri3_chase', 'vgfmri3_helper', 'vgfmri3_lemmings', 'vgfmri3_plaqueAttack', 'vgfmri3_zelda', ...
     'vgfmri3_chase + vgfmri3_helper + vgfmri3_lemmings + vgfmri3_plaqueAttack - vgfmri3_bait - vgfmri3_zelda', ... % agents - no agents
     'vgfmri3_lemmings + vgfmri3_plaqueAttack + vgfmri3_zelda - vgfmri3_bait - vgfmri3_chase - vgfmri3_helper'}, ... % shoot vs no shoot
     subjs);
