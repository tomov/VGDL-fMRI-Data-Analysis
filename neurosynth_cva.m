function neurosynth_cva(glmodel, use_smooth, lateralized, parcel_idx, cons, alpha)
    % c/p from neurosynth_rsa

%glmodel = 1;
%use_smooth = true;
%lateralized = true;
%parcel_idx = [1 2];
%cons = {'vgfmri3_bait', 'vgfmri3_chase', 'vgfmri3_helper', 'vgfmri3_lemmings', 'vgfmri3_plaqueAttack', 'vgfmri3_zelda'};
%alpha = 0.05;

% CVA with ROIs from neurosynth parcellation

printcode;

%rsa_idx = 1;
%lateralized = true;
%use_smooth = true;
%nperms = 1000; % 0 for no permutation tests

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:length(EXPT.subject);

% get ROIs
[roi_masks, region] = get_neurosynth_rois(lateralized);

if ~exist('parcel_idx', 'var') || isempty(parcel_idx)
    which = logical(ones(size(region)));

    tmp2 = join(cons, '+');
    filename = sprintf('mat/neurosynth_cva_%d_us=%d_l=%d_cons=%s_nroi=%d.mat', glmodel, use_smooth, lateralized, tmp2{1}, length(roi_masks));
else
    which = ismember(region, parcel_idx);
    roi_masks = roi_masks(which);

    tmp = join(cellfun(@num2str, num2cell(parcel_idx), 'UniformOutput', false), ',');
    tmp2 = join(cons, '+');
    filename = sprintf('mat/neurosynth_cva_%d_us=%d_l=%d_cons=%s_nroi=%d_pi=%s.mat', glmodel, use_smooth, lateralized, tmp2{1}, length(roi_masks), tmp{1});
end
filename

for i = 1:length(roi_masks)
    roi_mask = roi_masks{i};

    fprintf('ROI %d\n', i);

    tic
    CVA = ccnl_cva(EXPT, glmodel, cons, roi_mask, subjects);
    toc

    for s = 1:length(subjects)
        subj = subjects(s);

        CVAs(i,s).p = CVA{s}.p;
        CVAs(i,s).chi = CVA{s}.chi;
        CVAs(i,s).df = CVA{s}.df;
        CVAs(i,s).r = CVA{s}.r;
        CVAs(i,s).cva = CVA{s}.cva;
        CVAs(i,s).bic = CVA{s}.bic;
        CVAs(i,s).aic = CVA{s}.aic;

        ps(i,s,:) = CVA{s}.p;
    end

    if mod(i,10) == 0
        save(filename, '-v7.3');
    end
end

cnt = sum(ps < alpha, 3);
m = mean(cnt, 2);

filename
save(filename, '-v7.3');
