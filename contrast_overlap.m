% computes e part overlap of multiple contrasts

clear all;

glmodels = [85 51 52];
%glmodels = [86 82 83];
contrasts = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

EXPT = vgdl_expt();

% group-level settings
p = 0.01;
alpha = 0.05;
Dis = 20;
if ~exist('extent', 'var')
    extent = []; % use default cluster size
end
if ~exist('Num', 'var')
    Num = 1; % # peak voxels per cluster; default in bspmview is 3
end
direct = '+';
clusterFWEcorrect = true;

for i = 1:length(glmodels)
    glmodel = glmodels(i);
    contrast = contrasts{i};
    extent = [];

    [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent);

    if ~exist('cumulative_mask', 'var')
        cumulative_mask = C > 0;
    else
        cumulative_mask = cumulative_mask + (C>0);
    end
end

%ccnl_view_mask(cumulative_mask);
bspmview_wrapper(EXPT, cumulative_mask);

%
% calculate overlap for another contrast
%


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
