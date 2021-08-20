% get activations around given event(s)

EXPT = vgdl_expt();

subj_id = 1;
SPM_run_id = 1;

run_id = get_behavioral_run_id(subj_id, SPM_run_id);

dTRs = -2:6; % TRs relative to the event onset

% spherical mask around top ROI from contrast
glmodel = 21;
contrast = 'theory_change_flag';
mask_filenames = get_masks(glmodel, contrast, true, [], 3);
mask_filename = mask_filenames{1}; % top ROI only
[mask_format, mask, Vmask] = get_mask_format_helper(mask_filename);

% which events to extract time courses for
regs_fields = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};
visuals_fields = {'effects', 'avatar_collision_flag', 'new_sprites'};
onoff_fields = {'play_start', 'play_end'};
fields = [regs_fields, visuals_fields, onoff_fields];

% in lieu of get_regressors, get_visuals, get_onoffs, etc.
% because there's no mongo on NCF
load(sprintf('mat/get_regressors_subj%d_run%d.mat', subj_id, run_id), 'regs');
load(sprintf('mat/get_visuals_subj%d_run%d.mat', subj_id, run_id), 'visuals');
load(sprintf('mat/get_onoff_subj%d_run%d.mat', subj_id, run_id), 'onoff');

% event onsets, by type
onsets = struct;

% when we actually store the BOLD time courses for each event
activations = struct;

% get BOLD time course from ROI for given run
[Y, K, W, R, Y_run_id] = load_BOLD(EXPT, glmodel, subj_id, mask, Vmask);
Y_run = nanmean(Y(Y_run_id == SPM_run_id, :), 2);
assert(all(size(Y_run) == [EXPT.nTRs, 1]));

% extract event onsets
for i = 1:length(regs_fields)
    field = regs_fields{i};
    onsets.(field) = regs.timestamps(logical(regs.(field)));
end
for i = 1:length(visuals_fields)
    field = visuals_fields{i};
    onsets.(field) = visuals.timestamps(logical(visuals.(field)));
end
for i = 1:length(onoff_fields)
    field = onoff_fields{i};
    onsets.(field) = onoff.(field);
end

for i = 1:length(fields)
    field = fields{i};
    if ~isfield(activations, field)
        activations.(field) = [];
    end

    for j = 1:length(onsets.(field))
        TRs = round(onsets.(field)(j) / EXPT.TR) + dTRs;
        activations.(field) = [activations.(field); Y_run(TRs)'];
    end
end

save('PETHs.mat');

