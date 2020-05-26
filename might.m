% Run searchmight on the entire brain
%
% function might(rsa_idx, r, use_smooth, method, zsc)

clear all;
addpath(genpath('SearchmightToolbox.Darwin_i386.0.2.5/')); % for the cluster

rsa_idx = 1;
r = 4 / 1.5; % 4 mm
use_smooth = false;
method = 'gnb_searchmight'; % lda_ridge, lda_shrinkage, qda_shrinkage, svm_linear, svm_quadratic, svm_sigmoid, svm_rbf
zsc = 'none';


nruns = 6;

dirname = 'might';

if use_smooth
    EXPT = vgdl_expt();
    maskfile = 'masks/mask.nii';
    metafile = 'mat/might_meta.mat';
else
    EXPT = vgdl_expt_nosmooth();
    maskfile = 'masks/mask_nosmooth.nii';
    metafile = 'mat/might_meta_nosmooth.mat';
end

subjects = 1:length(EXPT.subject);

%use_tmaps = false; % <-- slightly better if true; but stick with betas for consistency w/ RDMs

%z_score = 'z-run';  % <-- better
%z_score = 'z-run-voxel';  % <-- nothing, though Storey's pFDR fucks up and gives q = 0.048 when none of them are actually significant
%z_score = 'z-none';   % <-- actually not bad
%z_score = 'z-manual'; % <-- hack; means manually z-score here; also nothing shows up
%z_score = 'z-random'; % <-- hack; for control, we use random activations

%classifier = 'gnb_searchmight'; % fast GNB; quick-n-dirty; gives weird axial dashes in accuracy maps, probably b/c of the way the scanner takes the images
%classifier = 'lda_shrinkage'; % based on Pereira and Botvinick (2010)


[mask] = ccnl_load_mask(maskfile);

%[meta] = createMetaFromMask(mask); <--- takes forever; load it instead (precomputed)
load(metafile);



dimx = size(mask, 1);
dimy = size(mask, 2); 
dimz = size(mask, 3); 
nvoxels = sum(mask(:));



% fix neighbors according to spherical searchlight
%
disp('computing neighbors');
tic
[meta.voxelsToNeighbours, meta.numberOfNeighbours] = might_computeNeighborsSphere(meta.colToCoord, r, mask);
toc




ams = nan(numel(subjects), nvoxels); % accuracy map for each subject 

alpha = 0.05; % for pFDR q-values
howmany = zeros(nvoxels, 1); % for each voxel, how many subjects have it significant (according to pFDR)

for s = 1:length(subjects)
    subj = subjects(s);
    fprintf('    subj %d\n', subj);


    disp('loading neural...');

    rsa = EXPT.create_rsa(rsa_idx, subj);

    % notice we load from nii every time, unlike RSA where we cache them
    % couple reasons:
    % 1) here use_smooth matters for which mask we use => need multiple caches
    % 2) this is not the bottleneck (as we only do this once)
    tic
    disp('loading betas from .nii files...');

    if rsa.use_beta_series
        % load betas
        %
        [B, names] = ccnl_get_beta_series(EXPT, rsa.glmodel, subj, rsa.event, mask);
        clear runs;
        for r = 1:nruns
            runs(contains(names, ['Sn(', num2str(r), ')']),:) = r;
        end
    else
        % load BOLD
        %
        [B, runs] = ccnl_get_activations(EXPT, rsa.glmodel, rsa.mask, subj, true, true); % whiten & filter; see Diedrichsen et al. 2016
        B = B{1};
        runs = runs{1};
    end
    toc

    % z score
    %

    switch zsc
        case 'time'
            B = zscore(B, 0, 1);

        case 'voxels'
            B = zscore(B, 0, 2);

        case 'run'
            for r = 1:nruns
                B(runs == r,:) = zscore(B(runs == r,:), 0, 1);
            end

        otherwise
            assert(isequal(zsc, 'none'));
    end


    % run searchlight
    %

    inputs = B;
    labels = rsa.model(1).features;
    labelsGroup = rsa.model(1).runs;

    disp('searchlighting...');

    tic
    [am,pm] = computeInformationMap(inputs,labels,labelsGroup,method,'searchlight', ...
                                    meta.voxelsToNeighbours,meta.numberOfNeighbours);
    ams(s,:) = am;
    toc


    % bookkeeping
    %

    [~, qm] = mafdr(pm'); % Storey (2002)
    howmany = howmany + (qm < alpha);

    disp('saving ouput');
    filename = sprintf('%s_accuracy_rsa=%d_subj=%d_folds=%d_r=%.4f_z%s_use_smooth=%d.nii', method, rsa_idx, subj, max(labelsGroup), r, zsc, use_smooth);
    % initialize an empty accuracy map
    [~, V, amap] = ccnl_load_mask(fullfile('masks', 'spmT_0001.nii'));
    V.fname = fullfile(dirname, filename); % change immediately!
    amap(:) = NaN; % clear

    % write accuracy map
    amap(mask) = am * 100;
    spm_write_vol(V, amap);

    % write thresholded map
    am(qm >= alpha) = NaN;
    amap(mask) = am * 100;
    V.fname = strrep(V.fname, '.nii', sprintf('_alpha=%.3f.nii', alpha));
    spm_write_vol(V, amap);

    % write p-value map
    filename = sprintf('%s_p-value_rsa=%d_subj=%d_folds=%d_r=%.4f_z%s_use_smooth=%d.nii', method, rsa_idx, subj, max(labelsGroup), r, zsc, use_smooth);
    V.fname = fullfile(dirname, filename);
    pmap = nan(size(amap));
    pmap(mask) = pm;
    spm_write_vol(V, pmap);

    % visualize
    %struc = fullfile('masks','mean.nii');
    %bspmview(V.fname, struc);
end


%% write map w/ # subjects for which voxel is significant (with pFDR) based on Pereira & Botvinick 2010
%

% initialize empty countmap 
filename = sprintf('%s_accuracy_countmap_%s_folds=%d_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii', classifier, event, max(labelsGroup), r, z_score, use_nosmooth, use_tmaps);
[~, V, countmap] = ccnl_load_mask(fullfile('masks', 'spmT_0001.nii'));
countmap(:) = NaN; % clear
V.fname = fullfile(dirname, filename); % change immediately!

% write countmap
%
[h, p, ci, stats] = ttest(ams, 1/6);
countmap(mask) = howmany;
spm_write_vol(V, countmap);

% visualize countmap
struc = fullfile('masks','mean.nii');
bspmview(V.fname, struc);

disp(V.fname);



%% write t-map based on Kriegeskorte & Bandettini 2007 
% note those have lots of negative t-values b/c "chance" is not 1/3 according to the classifier but slightly below it
% so don't use em
%

%{
% initialize empty tmap
filename = sprintf('%s_accuracy_tmap_%s_folds=%d_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii', classifier, event, max(labelsGroup), r, z_score, use_nosmooth, use_tmaps);
[~, V, tmap] = ccnl_load_mask(fullfile('masks', 'spmT_0001.nii'));
tmap(:) = NaN; % clear
V.fname = fullfile(dirname, filename); % change immediately!

% write tmap
% for each voxel, t-test subject accuracies against chance 
%
[h, p, ci, stats] = ttest(ams, 1/3);
tmap(mask) = stats.tstat;
spm_write_vol(V, tmap);

% visualize tmap
struc = fullfile('masks','mean.nii');
bspmview(V.fname, struc);

%}
