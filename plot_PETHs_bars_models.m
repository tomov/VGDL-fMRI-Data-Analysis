% plot results of gen_PETHs.m

close all;
clear all;


model_filenames = { ...
    fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP.mat'), ...
    fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_DQN.mat'), ...
};
model_names = {'EMPA', 'DQN'};
assert(length(model_filenames) == length(model_names));

for i = 1:length(model_filenames)
    i
    model_filenames{i}
    load(model_filenames{i}, 'mask_filenames', 'mask_name', 'regions', 'fields', 'activations', 'what', 'PETH_dTRs', 'EXPT');
    model_data(i).mask_filenames = mask_filenames;
    model_data(i).mask_name = mask_name;
    model_data(i).regions = regions;
    model_data(i).fields = fields;
    model_data(i).activations = activations;
    model_data(i).what = what;
    model_data(i).PETH_dTRs = PETH_dTRs;
    assert(isequal(model_data(1).mask_filenames, model_data(i).mask_filenames));
    assert(isequal(model_data(1).mask_name, model_data(i).mask_name));
    assert(isequal(model_data(1).regions, model_data(i).regions));
    assert(isequal(model_data(1).fields, model_data(i).fields));
end


field = 'theory_change_flag';

subjs = 1:1:32;

nROIs = length(mask_filenames);
nmodels = length(model_data);
nsubjects = length(subjs);

as = nan(nROIs,nmodels,nsubjects);

cmap = colormap(jet(nmodels));
t = PETH_dTRs * EXPT.TR; % s

% loop over masks
for m = 1:nROIs
    disp(mask_name{m});

    for i = 1:nmodels

        D = model_data(i).activations(m).(field)(subjs,:); % subj x TRs PETH's
        %as(m,i,:) = mean(D(:, PETH_dTRs > 5), 2); % average across time, ignoring baseline
        as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
    end
end

% Piggyback off of plot_gp_CV_rois.m
figure('pos', [49 329 2143 610]);

h = plot_gp_CV_rois_helper(as(:,:,:), 'ttest', 'mean', model_names, regions, 0, cmap, 5, 1:1);
if exist('what', 'var') && strcmp(what, 'GP')
    %title('Average Fisher z-transformed Pearson correlation change in ROIs');
    %ylabel('\Delta z');
    title('Average Fisher z-transformed Pearson correlation in ROIs');
    ylabel('z');
else
    title('Average BOLD change in ROIs');
    ylabel('\Delta BOLD');
end

% Prettyfy it 
% specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
%text(1.5, 0.75, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
%plot([2.5 2.5], [0 0.2], '--', 'color', [0.5 0.5 0.5]);
%text(3.5, 0.75, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
%plot([4.5 4.5], [0 0.2], '--', 'color', [0.5 0.5 0.5]);
%text(6, 0.75, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
%plot([7.5 7.5], [0 0.2], '--', 'color', [0.5 0.5 0.5]);
%text(8.5, 0.75, 'Early visual', 'fontsize', 12, 'HorizontalAlignment', 'center');
legend(fields(ix), 'interpreter', 'none');
