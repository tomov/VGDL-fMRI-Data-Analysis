
subj = 1;
mask = 'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii';

lateralized = true;
[roi_masks, region] = get_neurosynth_rois(lateralized);
mask = roi_masks{region == 38};

dist = {'cosine','correlation', 'euclidean'}

EXPTs = {vgdl_expt(), vgdl_expt_nosmooth()}
smooth = {'smooth', 'non-smooth'}


for sm = 1:2
    EXPT = EXPTs{sm};

    figure;
    idx = 1;

    A = ccnl_get_activations(EXPT, 1, mask, subj);

    for d = 1:3

        subplot(3,3,idx);
        idx = idx + 1;
        imagesc(squareRDMs(pdist(A, dist{d})));
        title('no z score');
        ylabel(dist{d});

        subplot(3,3,idx);
        idx = idx + 1;
        Y = zscore(A, 0, 1); % z score across time (like in ephys)
        imagesc(squareRDMs(pdist(Y, dist{d})));
        title('z score across time');

        subplot(3,3,idx);
        idx = idx + 1;
        Y = zscore(A, 0, 2); % z score across voxels (like in Baldassano et al. 2017)
        imagesc(squareRDMs(pdist(Y, dist{d})));
        title('z score across voxels');
    end

    figure;

    idx = 1;
    subplot(1,3,idx);
    idx = idx + 1;
    imagesc(A);
    title('no z score');
    ylabel('betas');

    subplot(1,3,idx);
    idx = idx + 1;
    Y = zscore(A, 0, 1); % z score across time (like in ephys)
    imagesc(Y);
    title('z score across time');

    subplot(1,3,idx);
    idx = idx + 1;
    Y = zscore(A, 0, 2); % z score across voxels (like in Baldassano et al. 2017)
    imagesc(Y);
    title('z score across voxels');
end

