% theory id sequence & unique theory HRRs => HRRs => Xx (convolved HRRs) => kernel
%
function [theory_kernel, theory_kernel_std, HRRs, Xx, r_id] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id, SPM, subsample_only)

    if ~exist('subsample_only', 'var')
        subsample_only = false;
    end

    %load('mat/SPM73.mat'); TOO SLOW!!!

    sigma_w = 1; % TODO param

    tot_1 = 0;
    tot_2 = 0;
    tot_3 = 0;

    nsamples = size(unique_theory_HRRs, 1); 

    clear Ks;
    for j = 1:nsamples

        unique_HRRs = squeeze(unique_theory_HRRs(j,:,:));

        % populate HRRs from unique_HRRs
        HRRs = nan(length(theory_id_seq), size(unique_HRRs,2));
        for i = 1:size(unique_HRRs,1)
            theory_id = i - 1; % b/c python is 0-indexed
            ix = find(theory_id_seq == theory_id);
            %HRRs(theory_id_seq == theory_id, :) = unique_HRRs(i,:); -- don't work
            for k = 1:length(ix)
                HRRs(ix(k), :) = unique_HRRs(i,:);
            end
        end

        %{
        % for sanity -- do with theory_change_flag, compare with GLM 3 
        for i = 1:size(HRRs,2)
            HRRs(:,i) = theory_change_flags;
        end
        %}

        if subsample_only
            [Xx, r_id] = subsample_HRRs(HRRs, ts, run_id, SPM);
        else
            [Xx, r_id] = convolve_HRRs(HRRs, ts, run_id, SPM);
        end

        Sigma_w = eye(size(Xx,2)) * sigma_w; % Sigma_p in Rasmussen, Eq. 2.4

        K = Xx * Sigma_w * Xx';

        if ~exist('Ks', 'var')
            Ks = nan(nsamples, size(K,1), size(K,2));
        end
        Ks(j,:,:) = K;
    end

    theory_kernel = squeeze(mean(Ks,1));
    theory_kernel_std = squeeze(std(Ks,0,1));
end

