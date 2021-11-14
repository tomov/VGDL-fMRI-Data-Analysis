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

load(fullfile(get_mat_dir(false), 'PETHs_atlas=AAL2_ungrouped2_BOLD.mat')); % !!!!!!!!!!!!
%ROI_ix = 1:length(mask_filenames);
ROI_ix = [      1      2      7     10     11     12     13     14     15]; 

mask_filenames = mask_filenames(ROI_ix);
mask_name = mask_name(ROI_ix);
regions = regions(ROI_ix);
activations = activations(ROI_ix);

figure('pos', [64 421 2282 838]);

% optionally plot theory change flag only
fields(find(strcmp(fields, 'theory_change_flag'))) = [];
%fields(find(strcmp(fields, 'sprite_change_flag'))) = [];
%fields(find(strcmp(fields, 'interaction_change_flag'))) = [];
%fields(find(strcmp(fields, 'termination_change_flag'))) = [];
fields(find(strcmp(fields, 'block_start'))) = [];
fields(find(strcmp(fields, 'block_end'))) = [];
fields(find(strcmp(fields, 'instance_start'))) = [];
fields(find(strcmp(fields, 'instance_end'))) = [];

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

    subplot(3, 3, m);
    hold on;

    for i = 1:nregressors
        field = fields{i};
        disp(field)

        D = activations(m).(field)(subjs,:); % subj x TRs PETH's
        as(m,i,:) = mean(D(:, PETH_dTRs > 0), 2); % average across time, ignoring baseline
    end
end

% Piggyback off of plot_gp_CV_rois.m
figure('position', [1147 521 1045 418]);
ix = 1:nregressors;
h = plot_gp_CV_rois_helper(as(:,ix,:), 'ttest', 'mean', fields(ix), regions, 0, cmap);
if exist('what', 'var') && strcmp(what, 'GP')
    title('Average Fisher z-transformed Pearson correlation change in ROIs');
    ylabel('\Delta z');
else
    title('Average BOLD change in ROIs');
    ylabel('\Delta BOLD');
end
