% plot results of gen_PETHs.m

close all;
clear all;

%load('mat/PETHs_glm=21_con=theory_change_flag_Num=3_sphere=4.0mm.mat')
%load('mat/PETHs_tag=tomov2018KL_sphere=10.0mm.mat')
%load('mat/PETHs_tag=hayley2021psi_sphere=10.0mm.mat');

%load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=10.0mm_.mat');
%load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=10.0mm__.mat');

%load('mat/PETHs_glm=21_con=theory_change_flag_odd_Num=1_sphere=6.0mm_.mat');

%load(fullfile(get_mat_dir(false), 'PETHs_glm=102_con=theory_change_flag_odd_Num=1_sphere=10.0mm.mat'));
%load(fullfile(get_mat_dir(false), 'PETHs_glm=102_con=theory_change_flag_odd_Num=1_sphere=10.0mm_GP.mat'));
%load(fullfile(get_mat_dir(true), 'PETHs_glm=102_con=theory_change_flag_odd_Num=1_sphere=10.0mm_.mat'));
% subselect ROIs
%ROI_ix = [1     2     3     5     7    11];

%load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_sphere=0.0mm_BOLD.mat'));
%load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL3v1_GP.mat'));

%load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GLM_102_BOLD.mat')); % !!!!!!!!!!!!
%ROI_ix = [      1      2      7     10     11     12     13     14     15]; 
load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP.mat')); % !!!!!!!!!!!!
%load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_GP_EMPA_GLM_102_GP_DQN.mat')); 
ROI_ix = 1:length(mask_filenames);

mask_filenames = mask_filenames(ROI_ix);
mask_name = mask_name(ROI_ix);
regions = regions(ROI_ix);
activations = activations(ROI_ix);

% optionally plot theory change flag only
%fields(find(strcmp(fields, 'theory_change_flag'))) = [];
fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
fields(find(strcmp(fields, 'termination_change_flag'))) = [];
fields(find(strcmp(fields, 'block_start'))) = [];
fields(find(strcmp(fields, 'block_end'))) = [];
fields(find(strcmp(fields, 'instance_start'))) = [];
fields(find(strcmp(fields, 'instance_end'))) = [];

%fields = {'theory_change_flag', 'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

subjs = 1:1:32;

nROIs = length(mask_filenames);
nregressors = length(fields);
nsubjects = length(subjs);

as = nan(nROIs,nregressors,nsubjects);

cmap = colormap(jet(length(fields)));
t = PETH_dTRs * EXPT.TR; % s

% loop over masks
for m = 1:nROIs
    disp(mask_name{m});

    for i = 1:nregressors
        field = fields{i};
        disp(field)

        D = activations(m).(field)(subjs,:); % subj x TRs PETH's
        as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
    end
end

% Piggyback off of plot_gp_CV_rois.m
figure('pos', [49 329 2143 610]);

ix = 1:nregressors;
h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap, 5, 1:3);
if exist('what', 'var') && strcmp(what, 'GP')
    title('Average Fisher z-transformed Pearson correlation change in ROIs');
    ylabel('\Delta z');
else
    title('Average BOLD change in ROIs');
    ylabel('\Delta BOLD');
end

% Prettyfy it 
% specifically for agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_GP_EMPA.mat');
text(1.5, 0.75, 'Frontal/Motor', 'fontsize', 12, 'HorizontalAlignment', 'center');
plot([2.5 2.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
text(3.5, 0.75, 'Dorsal/Parietal', 'fontsize', 12, 'HorizontalAlignment', 'center');
plot([4.5 4.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
text(6, 0.75, 'Ventral/Temporal', 'fontsize', 12, 'HorizontalAlignment', 'center');
plot([7.5 7.5], [0 0.8], '--', 'color', [0.5 0.5 0.5]);
text(8.5, 0.75, 'Early visual', 'fontsize', 12, 'HorizontalAlignment', 'center');
legend(fields(ix), 'interpreter', 'none');
