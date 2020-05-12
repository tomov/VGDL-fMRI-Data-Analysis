% export data for event segmentation in python BrainIAK

subj = 1;

mask = 'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii';
EXPT = vgdl_expt_nosmooth();
group_mask = 'masks/mask.nii';

nruns = 6;

filename = sprintf('mat/BOLD_for_eventseg_subj%d.mat', subj);

A = ccnl_get_activations(EXPT, 1, mask, subj);
BOLD = A{1};
assert(size(BOLD,1) == EXPT.nTRs * nruns);

[mask, V] = ccnl_load_mask(mask);

[x,y,z] = ind2sub(size(mask), find(mask));
coords = cor2mni([x y z], V.mat);

bounds = [];
HRF_lag = 6; % WTF TODO ask Chris Baldassano how they account for that 
for r = 1:nruns
    multi = vgdl_create_multi(3, subj, r);
    onsets = multi.onsets{1} + HRF_lag + (r - 1) * EXPT.run_duration;
    onsTRs = onsets / EXPT.TR;
    bounds = [bounds onsTRs];
end

save(filename, '-v6');
