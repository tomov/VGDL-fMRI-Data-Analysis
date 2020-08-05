function decode_gp_CV(subj_id)

%{
clear all;
close all;
subj_id = 1;
%}

use_smooth = true;
if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end
glmodel = 9;
maskfile = 'masks/ROI_x=42_y=28_z=26_1voxels_Sphere1.nii';


filename = sprintf('mat/decode_gp_CV_HRR_subj=%d.mat', subj_id);
filename

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
game_ids = [];
block_ids = [];

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

    game_ids = [game_ids; regs.game_ids];
    block_ids = [block_ids; regs.block_ids];
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

new_block_flag = block_ids(2:end) ~= block_ids(1:end-1);
new_block_flag = [0; new_block_flag];

% frames between which the theory is const
%
%bounds = union(find(interaction_change_flags), find(termination_change_flags));
bounds = find(theory_change_flags); % to be consistent
bounds = union(bounds, find(avatar_collision_flags));
bounds = union(bounds, find(new_block_flag));
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


load('mat/unique_HRR_subject_subj=1_K=10_N=10_E=0.050_nsamples=100_norm=1.mat', 'theory_HRRs', 'run_id', 'ts', 'theory_id_seq');
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

% init best theory sequence
%
n_theories = max(theory_id_seq);
r_best = r_orig;
nlz_best = nlz_orig;
theory_id_seq_best = theory_id_seq_orig;
nlz = nlz_orig; % for MCMC
theory_id_seq = theory_id_seq_orig; % for MCMC

% init theory id candidates for each game
%
tid_candidates = cell(max(game_ids),1);
for g = 1:max(game_ids)
    tid_candidates{g} = unique(theory_id_seq(game_ids == g));
end

% decode
%
%while (true) % repeat until we keep improving
for iter = 1:10000

    % pick random interval
    b = randi(length(bounds)-1);
    %for b = 1:length(bounds) - 1 % for each time point / set of consecutive time points
    st = bounds(b);
    en = bounds(b+1)-1;

    fprintf('between %d and %d\n', st, en);

    % pick theory id candidate at random
    assert(game_ids(st) == game_ids(en));
    assert(block_ids(st) == block_ids(en));
    g = game_ids(st);
    tid = randsample(tid_candidates{g}, 1);
    %for tid = 0:n_theories-1 % try assigning each candidate theory to that time point; note b/c of python, tid's are 0-indexed

    fprintf('     assign %d\n', tid);

    % candidate theory sequence = curreth theory, with one interval changed
    theory_id_seq_cand = theory_id_seq;
    theory_id_seq_cand(st:en) = tid;

    % get kernel for candidate theory sequence
    tic
    [theory_kernel, ~, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id_frames);
    toc

    tic
    ker = R*K*W*theory_kernel*W'*K'*R';
    toc

    % fit BOLD using CV and compare using r
    %
    %{
    tic
    [r_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, Y, ker, run_id_TRs);
    r = mean(r_CV, 1);
    toc

    % assess fit
    if r > r_best
        fprintf('             it''s better!! %d vs. %d\n', r, r_best);
        r_best = r;
        theory_id_seq_best = theory_id_seq;
    else
        theory_id_seq = theory_id_seq_best;
    end
    %}

    % fit BOLD w/o CV and compare using log marg lik
    %
    tic
    hyp = hyp_from_ker(nearestSPD(ker));
    hyp.lik = log(1); % TODO sigma_n = 1 const
    nlz_cand = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);
    toc

    % metropolis rule b/c proposal for tid is symmetric
    lme_cand = -nlz_cand;
    lme = -nlz;
    log_alpha = lme_cand - lme; % negative log marg liks
    u = rand();
    if log(u) <= log_alpha
        fprintf('             accept: %f vs. %f (log alpha = %f; %f <= %f\n', lme_cand, lme, log_alpha, u, exp(log_alpha));
        nlz = nlz_cand;
        theory_id_seq = theory_id_seq_cand;
    end

    % assess fit
    if nlz < nlz_best
        fprintf('             new best: %d vs. %d\n', nlz, nlz_best);
        nlz_best = nlz;
        theory_id_seq_best = theory_id_seq;
    end

    % optionally save
    if mod(iter,100) == 0 
        disp('saving');
        tic
        iter
        save(filename, '-v7.3');
        toc
    end
end

        
save(filename, '-v7.3');
disp('Done');


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


%{
figure;
imagesc(Xx);

figure;
imagesc(theory_kernel);
%}


