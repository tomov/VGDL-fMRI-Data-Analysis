function neurosynth_rsa(rsa_idx, use_smooth, lateralized, nperms, roi_idx_min, roi_idx_max, subbatch_size)

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

roi_idx_min = max(roi_idx_min, 1);
roi_idx_max = min(roi_idx_max, length(roi_masks));
roi_masks = roi_masks(roi_idx_min:roi_idx_max);

filename = sprintf('mat/neurosynth_rsa_%d_us=%d_l=%d_nroi=%d_nperms=%d_ri=%d-%d.mat', rsa_idx, use_smooth, lateralized, length(roi_masks), nperms, roi_idx_min, roi_idx_max);
disp(filename);

tic

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
        fprintf('perm = %d\n', i);
        
        EXPT.create_rsa = @(rsa_idx, subj_id) create_rsa(rsa_idx, subj_id, seed); % overwrite create_rsa with one that randomly shuffles 
        % run permuted RSA in (sub)batches
        %
        for s = 1:subbatch_size:length(roi_masks)
            e = min(length(roi_masks), s + subbatch_size - 1);
            fprintf('        perm s:e = %d:%d\n', s, e);

            [perm_Rhos(s:e,:,i)] = ccnl_rsa(EXPT, rsa_idx, roi_masks(s:e));
        end
    end

    save(filename, '-v7.3');

    % p-value = P(Rho >= rho | null)
    pval = mean(perm_Rhos >= Rho, 3);

    % adjusted rho = subtract mean & divide by stdev (for plotting)
    Rho_adj = (Rho - mean(perm_Rhos, 3)) ./ std(perm_Rhos, 0, 3);

    save(filename, '-v7.3');
end

disp(filename);

toc

% view RSA results
%
ccnl_rsa_view(EXPT, rsa_idx, 1, T, all_subject_rhos, roi_masks);
