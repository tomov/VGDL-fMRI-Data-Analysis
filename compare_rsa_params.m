
subj = 1;
mask = 'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii';

lateralized = true;
[roi_masks, region] = get_neurosynth_rois(lateralized);
mask = roi_masks{region == 38};

dist = {'cosine','correlation', 'euclidean', 'seuclidean', 'mahalanobis'}

EXPTs = {vgdl_expt(), vgdl_expt_nosmooth()}
smooth = {'smooth', 'non-smooth'}


for sm = 1:2
    EXPT = EXPTs{sm};

    figure;
    idx = 1;

    [A,runs] = ccnl_get_activations(EXPT, 1, mask, subj);
    A = A{1};
    runs = runs{1};

    for d = 1:5

        subplot(5,4,idx);
        idx = idx + 1;
        imagesc(squareRDMs(pdist(A, dist{d})));
        title('no z score');
        ylabel(dist{d});

        subplot(5,4,idx);
        idx = idx + 1;
        Y = zscore(A, 0, 1); % z score across time (like in ephys)
        imagesc(squareRDMs(pdist(Y, dist{d})));
        title('z score across time');

        subplot(5,4,idx);
        idx = idx + 1;
        for r = 1:6
            Y(runs == r,:) = zscore(A(runs == r,:), 0, 1); % z score across time, within run (like in ephys)
        end
        imagesc(squareRDMs(pdist(Y, dist{d})));
        title('z score within run');

        if ~isequal(dist{d}, 'mahalanobis')
            subplot(5,4,idx);
            idx = idx + 1;
            Y = zscore(A, 0, 2); % z score across voxels (like in Baldassano et al. 2017)
            imagesc(squareRDMs(pdist(Y, dist{d})));
            title('z score across voxels');
        end
    end

    figure;

    idx = 1;
    subplot(1,4,idx);
    idx = idx + 1;
    imagesc(A);
    title('no z score');
    ylabel('betas');

    subplot(1,4,idx);
    idx = idx + 1;
    Y = zscore(A, 0, 1); % z score across time (like in ephys)
    imagesc(Y);
    title('z score across time');

    subplot(1,4,idx);
    idx = idx + 1;
    for r = 1:6
        Y(runs == r,:) = zscore(A(runs == r,:), 0, 1); % z score across time, within run (like in ephys)
    end
    imagesc(Y);
    title('z score within run');

    subplot(1,4,idx);
    idx = idx + 1;
    Y = zscore(A, 0, 2); % z score across voxels (like in Baldassano et al. 2017)
    imagesc(Y);
    title('z score across voxels');
end

