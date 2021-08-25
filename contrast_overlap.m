% computes e part overlap of multiple contrasts

clear all;

glmodels = [85 51 52];
contrasts = {'sprite_change_flag', 'interaction_change_flag', 'termination_change_flag'};

EXPT = vgdl_expt();

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
direct = '+/-';
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
