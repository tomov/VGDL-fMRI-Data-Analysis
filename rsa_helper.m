function rsa_helper(EXPT, rsa_idx, roi_masks, filename, nperms, subbatch_size, region, which)

    % do the actual RSA given binary masks


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

                % USE ALL SUBJECT RHOS !!! x8 data points!! 
                % JUST KIDDING DON'T DO IT!! would you test for significance for subject A by doing a shuffle on subject B? no -- b/c some subjects might have e.g. high activations overall, leading to more similar correlations overall (see Walther 2015)
                [perm_Rhos(s:e,:,i)] = ccnl_rsa(EXPT, rsa_idx, roi_masks(s:e));
            end

            if mod(i,20) == 0
                save(filename, '-v7.3');
            end
        end
                
        save(filename, '-v7.3');

        % p-value = P(Rho >= rho | null)
        pval = mean(perm_Rhos >= Rho, 3); % TODO include Rho as well

        % adjusted rho = subtract mean & divide by stdev (for plotting)
        Rho_adj = (Rho - mean(perm_Rhos, 3)) ./ std(perm_Rhos, 0, 3);

        save(filename, '-v7.3');
    end

    filename

    % view RSA results
    %
    ccnl_rsa_view(EXPT, rsa_idx, 1, T, all_subject_rhos, roi_masks);
