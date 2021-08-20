function [masks, region] = get_masks(glmodel, contrast, clusterFWEcorrect, extent, Num, sphere)

    % copied from Exploration

    if ~exist('sphere', 'var')
        sphere = 10; % 10 mm radius by default
    end

    EXPT = vgdl_expt();


    % it's an actual contrast
    %

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

    [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent);

    r = sphere / 1.5; % 10 mm radius

    % create spherical masks around peak voxel of each cluster (intersected with cluster)
    %
    for c = 1:length(region)
        masks{c} = sprintf('sphere_glm%d_%s_%d_%d_%d_r=%dmm.nii', glmodel, replace(contrast, ' ', '_'), mni(c,1), mni(c,2), mni(c,3), round(r * 1.5));
        %cmask = CI == CI(cor(c,1), cor(c,2), cor(c,3));
        ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r, masks{c}, []);
    end

