% computes e part overlap of multiple contrasts

clear all;

load(fullfile(get_mat_dir(), 'agg_gp_CV_us=1_glm=1_model=EMPA_theory_nsamples=100_project=1_fast=1.mat'));
tmap_filename = bspmview_save_map(EXPT, tmap);

EXPTs = {vgdl_expt(), tmap_filename}; % Take advantage of hack that allows us to plot any custom nmap
glmodels = [102 , -1];
contrasts = {'theory_change_flag', ''};

%glmodels = [85 51 52];
%glmodels = [86 82 83];
%contrasts = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};


% group-level settings
p = 0.001;
alpha = 0.05;
Dis = 20;
if ~exist('extent', 'var')
    extent = []; % use default cluster size
end
if ~exist('Num', 'var')
    Num = 1; % # peak voxels per cluster; default in bspmview is 3
end
direct = '+';
df = 31;
clusterFWEcorrect = true;

for i = 1:length(glmodels)
    EXPT = EXPTs{i};
    glmodel = glmodels(i);
    contrast = contrasts{i};
    extent = [];

    [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent, df);

    if ~exist('cumulative_mask', 'var')
        cumulative_mask = C > 0;
    else
        cumulative_mask = cumulative_mask + (C>0);
    end
end

%ccnl_view_mask(cumulative_mask);
bspmview_wrapper(vgdl_expt(), cumulative_mask);

%
% calculate overlap for another contrast
%


%{
p = 0.001;
glmodel = 21;
contrast = 'theory_change_flag';
extent = [];
clusterFWEcorrect = true;
[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent);

num_overlap = [];
for i = 1:size(cor, 1)
    num_overlap(i,:) = cumulative_mask(cor(i,1), cor(i,2), cor(i,3));
end

table(region, mni, stat, extent, num_overlap)
%}
