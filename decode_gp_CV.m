
clear all;
close all;

%load('../glmOutput/model3/subj1/multi1.mat');

subj_id = 1;

%
% from vgdl_create_multi.m
%

conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')

timestamps = [];
collision_flags = [];
avatar_collision_flags = [];
interaction_change_flags = [];
termination_change_flags = [];
theory_change_flags = [];

run_id_2 = [];
for r = 1:6
    query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, r) % in python we index runs from 0 (but not subjects) 

    run = find(conn, 'runs', 'query', query)
    assert(length(run) == 1);

    regs = get_regressors(subj_id, run, conn, true);

    [~, visuals] = get_visuals(subj_id, run, conn, true);
    assert(length(visuals.durations) == length(regs.timestamps));

    timestamps = [timestamps; regs.timestamps];

    tc = regs.theory_change_flag;
    tc = logical(tc);
    assert(all(regs.timestamps(tc) == regs.theory_change_flag_onsets));

    ef = visuals.effectsByCol;
    ef = ef > 0;
    collision_flags = [collision_flags; ef];

    av = visuals.avatar_collision_flag;
    av = logical(av);
    avatar_collision_flags = [avatar_collision_flags; av];

    theory_change_flags = [theory_change_flags; regs.theory_change_flag];
    interaction_change_flags = [interaction_change_flags; regs.interaction_change_flag];
    termination_change_flags = [termination_change_flags; regs.termination_change_flag];

    run_id_2 = [run_id_2; ones(length(regs.theory_change_flag),1) * r];
    % theory_change_flag is unreliable -- sprite types flicker randomly, 
    % so sanity check with other flags
    %assert(all(ismember(find(regs.interaction_change_flag), find(visuals.effectsByCol))));
    %assert(all(ismember(find(regs.termination_change_flag), find(visuals.effectsByCol))));
end
theory_change_flags = logical(theory_change_flags);
interaction_change_flags = logical(interaction_change_flags);
termination_change_flags = logical(termination_change_flags);


% frames between which the theory is const
%
bounds = union(find(interaction_change_flags), find(termination_change_flags))
bounds = union(bounds, find(avatar_collision_flags))
%bounds = union(bounds, find(collision_flags)) % -- too many...



% actually populate HRRs from unique HRRs (theory id -> HRR) and theory id sequence
%


load('mat/unique_HRR_subject_subj=1_K=10_N=10_E=0.050_nsamples=10_norm=1.mat');
run_id = run_id';
assert(all(run_id == run_id_2));
ts = ts';
assert(immse(timestamps, ts) < 1e-20);
assert(length(ts) == length(theory_id_seq));

% sanity
for r = 1:6
    multi = vgdl_create_multi(3, 1, r);
    assert(immse(ts(theory_change_flags & run_id == r), multi.onsets{1}) < 1e-20);

    load(sprintf('../glmOutput/model3/subj1/multi%d.mat', r));
    assert(immse(ts(theory_change_flags & run_id == r), onsets{1}) < 1e-20);
end

nsamples = 10;

sigma_w = 1; % TODO param

tot_1 = 0;
tot_2 = 0;
tot_3 = 0;

load('mat/SPM73.mat');

clear Ks;
for j = 1:nsamples

    j
    tic;

    unique_HRRs = squeeze(theory_HRRs(j,:,:));

    % populate HRRs from unique_HRRs
    HRRs = nan(length(ts), size(unique_HRRs,2));
    for i = 1:size(unique_HRRs,1)
        theory_id = i - 1; % b/c python is 0-indexed
        ix = find(theory_id_seq == theory_id);
        %HRRs(theory_id_seq == theory_id, :) = unique_HRRs(i,:); -- don't work
        for k = 1:length(ix)
            HRRs(ix(k), :) = unique_HRRs(i,:);
        end
    end

    tot_1 = tot_1 + toc;
    tic;

    %{
    % for sanity -- do with theory_change_flag, compare with GLM 3 
    for i = 1:size(HRRs,2)
        HRRs(:,i) = theory_change_flags;
    end
    %}

    [Xx, r_id] = convolve_HRRs(HRRs, ts, run_id, SPM);

    tot_2 = tot_2 + toc;
    tic;

    Sigma_w = eye(size(Xx,2)) * sigma_w; % Sigma_p in Rasmussen, Eq. 2.4

    K = Xx * Sigma_w * Xx';

    if ~exist('Ks', 'var')
        Ks = nan(nsamples, size(K,1), size(K,2));
    end
    Ks(j,:,:) = K;

    tot_3 = tot_3 + toc;
end

tot_1
tot_2
tot_3

theory_kernel = squeeze(mean(Ks,1));
theory_kernel_std = squeeze(std(Ks,0,1));




use_smooth = true;
glmodel = 9;
maskfile = 'masks/ROI_x=42_y=28_z=26_1voxels_Sphere1.nii';
[r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, maskfile, theory_kernel);

r_CV_1 = r_CV;
R2_CV_1 = R2_CV;
load('mat/fit_gp_CV_HRR_subj=1_us=1_glm=9_mask=mask_theory.mat', 'r_CV', 'R2_CV', 'mask');
whole_brain_mask = mask;
[mask] = ccnl_load_mask(maskfile);
which = mask(whole_brain_mask);
r_CV = r_CV(:,which);
R2_CV = R2_CV(:,which);

%{
% for sanity -- after convolning theory_change_flag manually, compare with GLM 3 

load('../glmOutput/model3/subj1/SPM.mat');

Xt = sum(SPM.xX.X(:,contains(SPM.xX.name, 'theory_change_flag')),2);


figure;
subplot(1,2,1);
imagesc(Xt);
title('from SPM');
subplot(1,2,2);
imagesc(Xx);
title('manual');
%}


figure;
imagesc(Xx);

figure;
imagesc(theory_kernel);



%
% from convolve_HRRs() in HRR.py
%
function [Xx, r_id] = convolve_HRRs(HRRs, ts, run_id, SPM)

    nruns = length(SPM.nscan);
    assert(nruns == 6);

    Xx = [];
    r_id = [];

    for s = 1:nruns
        
        which = run_id == s;

        % from spm_get_ons.m
        %
        k = SPM.nscan(s);
        assert(k == 283);

        T = SPM.xBF.T;
        assert(T == 16);

        dt = SPM.xBF.dt;
        assert(dt == 0.1250);

        UNITS = SPM.xBF.UNITS;
        assert(isequal(SPM.xBF.UNITS, 'secs'));
        TR = 1;

        bf = SPM.xBF.bf;
        assert(size(bf,1) == 257);

        % from spm_get_ons.m
        %

        ons = ts(which);
        u = HRRs(which,:);
        ton = round(ons*TR/dt) + 33; % 32 bin offset
        %toff = ton + 1; % for sanity, with impulse regressors only, e.g. theory_change_flag w/ GLM 3
        toff = ton(2:end); % frames are back-to-back
        toff = [toff; ton(end) + T]; % 1 s duration for last one, to match other between-block / level durations TODO more rigorous
        sf = zeros((k*T + 128), size(u,2));

        assert(all(ton >= 0));
        assert(all(ton < size(sf,1)));
        assert(all(toff >= 0));
        assert(all(toff < size(sf,1)));

        for j = 1:length(ton)
            sf(ton(j),:) = sf(ton(j),:) + u(j,:);
            sf(toff(j),:) = sf(toff(j),:) - u(j,:);
        end
        
        sf = cumsum(sf);
        sf = sf(1:(k*T + 32),:);                 %  32 bin offset

        % from spm_Volterra.m
        %

        %{
        % don't use convolution matrix; it's slower...
        if ~exist('A', 'var')
            A = convmtx(bf, size(sf,1));
            d = 1:size(sf,1);
        end
        X = A * sf;
        X = X(d,:);
        %}
        X = zeros(size(sf));
        for i = 1:size(sf,2)
            x = sf(:,i);
            d = 1:length(x);
            x = conv(x, bf);
            x = x(d);
            X(:,i) = x;
        end

        % from spm_fMRI_design.m
        %

        fMRI_T = SPM.xBF.T;
        fMRI_T0 = SPM.xBF.T0;

        assert(k == 283);
        idx = (0:(k - 1))*fMRI_T + fMRI_T0 + 32;
        X = X(idx,:);

        Xx = [Xx; X];
        r_id = [r_id; ones(size(X,1),1) * s];
    end

end


