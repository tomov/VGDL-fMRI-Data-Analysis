function neurosynth_rsa(rsa_idx, use_smooth, lateralized, nperms, parcel_idx, subbatch_size)

% copied from Exploration / neurosynth_CV.m

printcode;

%rsa_idx = 1;
%lateralized = true;
%use_smooth = true;
%nperms = 1000; % 0 for no permutation tests

if ~exist('subbatch_size', 'var')
    subbatch_size = 50; % don't do all ROIs at once; we OOM b/c all the Neural RDMs; too few though, and you end up loading betas too often
end

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

% get ROIs
[roi_masks, region] = get_neurosynth_rois(lateralized);

if ~exist('parcel_idx', 'var') || isempty(parcel_idx)
    filename = sprintf('mat/neurosynth_rsa_%d_us=%d_l=%d_nperms=%d_nroi=%d.mat', rsa_idx, use_smooth, lateralized, nperms, length(roi_masks));
else
    which = ismember(region, parcel_idx);
    roi_masks = roi_masks(which);

    tmp = join(cellfun(@num2str, num2cell(parcel_idx), 'UniformOutput', false), ',');
    filename = sprintf('mat/neurosynth_rsa_%d_us=%d_l=%d_nperms=%d_nroi=%d_pi=%s.mat', rsa_idx, use_smooth, lateralized, nperms, length(roi_masks), tmp{1});
end
filename

% run RSA in (sub)batches
%
for s = 1:subbatch_size:length(roi_masks)
    e = min(length(roi_masks), s + subbatch_size - 1);
    fprintf('s:e = %d:%d\n', s, e);

    [Rho(s:e,:), H(s:e,:), T(s:e,:), P(s:e,:), all_subject_rhos(s:e,:,:)] = ccnl_rsa(EXPT, rsa_idx, roi_masks(s:e));
end

Behavioral = ccnl_behavioral_rdms(EXPT, rsa_idx); % for plotting

save(filename, '-v7.3');

% permutations
%
if nperms > 0
    rng('shuffle');

    create_rsa = EXPT.create_rsa;

    perm_Rhos = nan(size(Rho,1), size(Rho,2), nperms);
    for i = 1:nperms
        seed = randi(1000000); % notice this will get affected by the rng calls inside create_rsa, but that's okay
        fprintf('perm = %d, seed = %d\n', i, seed);
        
        EXPT.create_rsa = @(rsa_idx, subj_id) create_rsa(rsa_idx, subj_id, seed); % overwrite create_rsa with one that randomly shuffles 
        % run permuted RSA in (sub)batches
        %
        for s = 1:subbatch_size:length(roi_masks)
            e = min(length(roi_masks), s + subbatch_size - 1);
            fprintf('        perm s:e = %d:%d\n', s, e);

            [perm_Rhos(s:e,:,i)] = ccnl_rsa(EXPT, rsa_idx, roi_masks(s:e));
        end

        if mod(i,100) == 0
            save(filename, '-v7.3');
        end
    end
            
    save(filename, '-v7.3');

    % p-value = P(Rho >= rho | null)
    pval = mean(perm_Rhos >= Rho, 3);

    % adjusted rho = subtract mean & divide by stdev (for plotting)
    Rho_adj = (Rho - mean(perm_Rhos, 3)) ./ std(perm_Rhos, 0, 3);

    save(filename, '-v7.3');
end

filename

% view RSA results
%
ccnl_rsa_view(EXPT, rsa_idx, 1, T, all_subject_rhos, roi_masks);
