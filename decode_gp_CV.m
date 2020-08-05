
clear all;
close all;

%load('../glmOutput/model3/subj1/multi1.mat');

subj_id = 1;

%
% from vgdl_create_multi.m
%

try
    conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')
catch
    disp('NCF....................');
end


%
% extract event sequence, to form bounds where we will keep the theory constant
%

timestamps = [];
collision_flags = [];
avatar_collision_flags = [];
interaction_change_flags = [];
termination_change_flags = [];
theory_change_flags = [];

run_id_2 = [];
for r = 1:6
    query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, r); % in python we index runs from 0 (but not subjects) 

    try
        run = find(conn, 'runs', 'query', query)
        assert(length(run) == 1);

        regs = get_regressors(subj_id, run, conn, true);

        [~, visuals] = get_visuals(subj_id, run, conn, true);

    catch e
        e
        fprintf('loading from cache...');
        load(sprintf('mat/get_regressors_subj%d_run%d.mat', subj_id, r));
        load(sprintf('mat/get_visuals_subj%d_run%d.mat', subj_id, r));
    end
        
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
bounds = union(find(interaction_change_flags), find(termination_change_flags));
bounds = union(bounds, find(avatar_collision_flags));
bounds = [1; bounds; length(theory_change_flags) + 1];
%bounds = union(bounds, find(collision_flags)) % -- too many...


%{
% optionally remove 1-frame short intervals by concatenating them
%
clear_which = logical(zeros(size(bounds)));
prev_len = nan;
for i = 1:length(bounds)-1
    len = bounds(i+1) - bounds(i);
    if len <= 4 && prev_len <= 4
        clear_which(i) = 1;
    end
    prev_len = len;
end
bounds(clear_which) = [];
%}


% actually populate HRRs from unique HRRs (theory id -> HRR) and theory id sequence
%


load('mat/unique_HRR_subject_subj=1_K=10_N=10_E=0.050_nsamples=10_norm=1.mat', 'theory_HRRs', 'run_id', 'ts', 'theory_id_seq');
unique_theory_HRRs = theory_HRRs;
run_id_frames = run_id';
assert(all(run_id_frames == run_id_2));
ts = ts';
assert(immse(timestamps, ts) < 1e-20);
assert(length(ts) == length(theory_id_seq));

% sanity
for r = 1:6
    multi = vgdl_create_multi(3, 1, r);
    assert(immse(ts(theory_change_flags & run_id_frames == r), multi.onsets{1}) < 1e-20);

    %{
    % for local sanity check
    load(sprintf('../glmOutput/model3/subj1/multi%d.mat', r));
    assert(immse(ts(theory_change_flags & run_id == r), onsets{1}) < 1e-20);
    %}
end


% create kernel from theory id sequence
%
[theory_kernel, ~, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id_frames);

% pre-load BOLD
%
use_smooth = true;
if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end
glmodel = 9;
maskfile = 'masks/ROI_x=42_y=28_z=26_1voxels_Sphere1.nii';
[mask_format, mask, Vmask] = get_mask_format_helper(maskfile);
[Y, K, W, R, run_id_TRs] = load_BOLD(EXPT, glmodel, subj_id, mask, Vmask);
Y = R*K*W*Y;
assert(size(Y,2) == 1); % single voxel
y = Y;

% init GP stuff
%
addpath(genpath('/ncf/gershman/Lab/scripts/gpml'));

n = size(y, 1);
x = [1:n]';
meanfun = @meanConst;
covfun = {@covDiscrete, n};
likfun = @likGauss;


% fit GP with original kernel
%
%[r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, maskfile, theory_kernel);
ker = R*K*W*theory_kernel*W'*K'*R';
[r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, Y, ker, run_id_TRs);

% fit manually too -- faster?
hyp = hyp_from_ker(nearestSPD(ker));
hyp.lik = log(1); % TODO sigma_n = 1 const
nlz = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);

assert(size(r_CV,2) == 1);
r_orig = mean(r_CV, 1);
nlz_orig = nlz;
theory_id_seq_orig = theory_id_seq;
nlz_orig
r_orig

% decode
%
n_theories = max(theory_id_seq);
r_best = r_orig;
nlz_best = nlz_orig;
theory_id_seq_best = theory_id_seq_orig;

while (true) % repeat until we keep improving

    %for b = 1:length(bounds) - 1 % for each time point / set of consecutive time points
    for b = length(bounds)-1:-1:1 % for each time point / set of consecutive time points
        st = bounds(b);
        en = bounds(b+1)-1;

        fprintf('between %d and %d\n', st, en);

        for tid = 1:n_theories % try assigning each candidate theory to that time point
            fprintf('     assign %d\n', tid);

            % change theory 
            theory_id_seq(st:en) = tid - 1; % python

            % get kernel
            tic
            [theory_kernel, ~, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id_frames);
            toc

            % fit BOLD
            %[r_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, maskfile, theory_kernel);
            tic
            ker = R*K*W*theory_kernel*W'*K'*R';
            toc

            %{
            tic
            [r_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, Y, ker, run_id_TRs);
            r = mean(r_CV, 1);
            toc
            %}

            tic
            clear r;
            hyp = hyp_from_ker(nearestSPD(ker));
            hyp.lik = log(1); % TODO sigma_n = 1 const
            nlz = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);
            toc

            % assess fit
            %if r > r_best
            if nlz < nlz_best
                %fprintf('             it''s better!! %d vs. %d\n', r, r_best);
                fprintf('             it''s better!! %d vs. %d\n', nlz, nlz_best);
                %r_best = r;
                nlz_best = nlz;
                theory_id_seq_best = theory_id_seq;
            end
        end

        theory_id_seq = theory_id_seq_best;
    end

end



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
% theory id sequence & unique theory HRRs => HRRs => Xx (convolved HRRs) => kernel
%
function [theory_kernel, theory_kernel_std, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id)

    load('mat/SPM73.mat');

    sigma_w = 1; % TODO param

    tot_1 = 0;
    tot_2 = 0;
    tot_3 = 0;

    nsamples = size(unique_theory_HRRs, 1); 

    clear Ks;
    for j = 1:nsamples

        tic;

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

    %tot_1
    %tot_2
    %tot_3

    theory_kernel = squeeze(mean(Ks,1));
    theory_kernel_std = squeeze(std(Ks,0,1));
end



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


