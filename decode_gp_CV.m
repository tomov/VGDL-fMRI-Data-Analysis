%function decode_gp_CV(subj_id)

clear all;
close all;
subj_id = 1;

use_smooth = true;
if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end
glmodel = 9;
maskfile = 'masks/ROI_x=42_y=28_z=26_1voxels_Sphere1.nii';


filename = sprintf('mat/decode_gp_CV_HRR_subj=%d_minint=300_2.mat', subj_id);
filename

%
% from vgdl_create_multi.m
%

try
    conn = mongo('holy7c22105.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54', 'UserName', 'reader', 'Password', 'parolatamadafaqa')
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
run_ids = [];
for r = 1:6  % TODO hardcoded
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
    run_ids = [run_ids; ones(length(regs.theory_change_flag),1) * r];

    theory_change_flags = [theory_change_flags; regs.theory_change_flag];
    interaction_change_flags = [interaction_change_flags; regs.interaction_change_flag];
    termination_change_flags = [termination_change_flags; regs.termination_change_flag];

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

bounds = cleanup_bounds(bounds);

%figure;
%subplot(2,1,1);
%plot(bounds_orig);
%subplot(2,1,2);
%plot(bounds);



% actually populate HRRs from unique HRRs (theory id -> HRR) and theory id sequence
%


load(sprintf('mat/unique_HRR_subject_subj=%d_K=10_N=10_E=0.050_nsamples=100_norm=1.mat', subj_id), 'theory_HRRs', 'run_id', 'ts', 'theory_id_seq', 'play_key', 'gameStrings', 'unique_theories_filename');
unique_theory_HRRs = theory_HRRs;
unique_theory_HRRs = unique_theory_HRRs(1:20,:,:); % TODO !!!!!!!!!!!!!!!!
run_id_frames = run_id';
assert(all(run_id_frames == run_ids));
ts = ts';
assert(immse(timestamps, ts) < 1e-20);
assert(length(ts) == length(theory_id_seq));

% sanity
for r = 1:6
    multi = vgdl_create_multi(3, subj_id, r);
    assert(immse(ts(theory_change_flags & run_id_frames == r), multi.onsets{1}) < 1e-20);

    % for local sanity check
    %load(sprintf('../glmOutput/model3/subj%d/multi%d.mat', subj_id, r));
    %assert(immse(ts(theory_change_flags & run_id == r), onsets{1}) < 1e-20);
end




tc = find(theory_id_seq(2:end) ~= theory_id_seq(1:end-1)) + 1; % frame indices of theory changes 
tc = [1 tc length(theory_id_seq)+1];
bounds = union(bounds, tc);

blk = find(block_ids(2:end) ~= block_ids(1:end-1)) + 1; % frame after block change
bounds = union(bounds, blk);


% create kernel from theory id sequence
%

load('mat/SPM73.mat');

tic
[theory_kernel, ~, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id_frames, SPM);
toc



% pre-load BOLD
%
[mask_format, mask, Vmask] = get_mask_format_helper(maskfile);
%[Y, K, W, R, run_id_TRs] = load_BOLD(EXPT, glmodel, subj_id, mask, Vmask);
load mat/load_BOLD_39.mat % TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Y = R*K*W*Y;
%assert(size(Y,2) == 1); % single voxel
y = Y; % TODO rename to Y

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
[r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, Y, ker, run_id_TRs, x, y, meanfun, covfun, likfun);

% fit manually too -- faster? for 1 voxel only
%hyp = hyp_from_ker(nearestSPD(ker));
%hyp.lik = log(1); % TODO sigma_n = 1 const
%nlz = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);
nlz = fit_gp_simple(Y, ker, x, y, meanfun, covfun, likfun);
nlz = mean(nlz);

%assert(size(r_CV,2) == 1);
r_CV = mean(r_CV, 2);
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

nlz = nlz_orig; % for greedy
theory_id_seq = theory_id_seq_orig; % for greedy

% init theory id candidates for each game
%
tid_candidates = cell(max(game_ids),1);
for g = 1:max(game_ids)
    tid_candidates{g} = unique(theory_id_seq(game_ids == g));
end

% decode
%



%
% step 1) first fit the theory change flag points
% for each change point, greedily try moving it left or right, until it stops improving fit
% pick the side that gives best fit
% rationale: GP is kind of like RSA; change points would be more salient (RDM blocks)
% even if the underlying theories are slightly incorrect
% also note that these changes are mostly local, as in they are mostly affected by (and affecting) what's going on nearby
%

figure;

t_next = nan;
t_prev = nan;
% go backwards, b/c 
% 1) we have more data (BOLD) for later intervals b/c they are larger b/c EMPA (and humans?) have more or less stopped learning
% 2) we are more certain about the later theories b/c of convergence / prior doesn't matter as much
for i = length(tc)-1:-1:2 %
%for i = 2:1:length(tc)-1 % not as good of a fit, and it makes sense to go in reverse direction

    %i = randi(length(tc)-2) + 1;
    i

    % NOTE: do this only if iterating in reverse
    if isnan(t_next)
        t_next = tc(i+1); % next theory change
    else
        t_next = t; % from last iteration
    end
    t = tc(i); % current theory_change frame; we will try moving it around
    %t_next = tc(i+1); % next theory change NOTE: do this only if iterating forward
    t_prev = tc(i-1); % previous theory change; we don't want to go past it, NOTE: do this only if iterating in reverse
    %{
    % NOTE: do this only if iterating forward
    if isnan(t_prev)
        t_prev = tc(i-1); % next theory change
    else
        t_prev = t; % from last iteration
    end
    %}

    % keep pushing change theory change point back, until we hit the previous change point
    % or the start of the block (so we don't assign the same theory across blocks)
    % note that if we hit the previous change point, we push it back on the next iteration of the outer loop
    %
    j = find(bounds == t);
    assert(length(j) == 1);

    fprintf('LEFT: b = %d, t_prev = %d, b_prev = %d\n', bounds(j), t_prev, bounds(j-1));

    nlz_left = nlz;
    j_left = j;
    theory_id_seq_left = theory_id_seq;
    theory_id_seq_cand = theory_id_seq;
    while j > 1 && bounds(j) > t_prev
        j = j - 1;
        if block_ids(bounds(j)) ~= block_ids(t)
            break
        end
        theory_id_seq_cand(bounds(j):t-1) = theory_id_seq(t);

        tic
        nlz_cand = get_theory_seq_nlz(unique_theory_HRRs, theory_id_seq_cand, ts, run_id_frames, R, K, W, x, y, meanfun, covfun, likfun, SPM);
        toc

        fprintf('%d -> %d: nlz = %.5f vs. %.5f\n', t, bounds(j), nlz_cand, nlz_left);
        if nlz_cand <= nlz_left
            fprintf('                        new best: %d vs. %d\n', nlz_cand, nlz_left);
            nlz_left = nlz_cand;
            j_left = j;
            theory_id_seq_left = theory_id_seq_cand;
        else
            disp('                                          not better');
            break;
        end

        fprintf('%d -> %d\n', t, bounds(j));
        plot(theory_id_seq);
        hold on;
        plot(block_ids * 100 + 100);
        plot(bounds(j), theory_id_seq(t), '*', 'color', 'red');
        hold off;
        drawnow
        %pause(1);
    end


    % same but go right
    %
    j = find(bounds == t);
    assert(length(j) == 1);

    fprintf('RIGHT: b = %d, t_next = %d, b_next = %d\n', bounds(j), t_next, bounds(j+1));

    nlz_right = nlz;
    j_right = j;
    theory_id_seq_right = theory_id_seq;
    theory_id_seq_cond = theory_id_seq;
    while j < length(bounds) && bounds(j) < t_next
        j = j + 1;
        if bounds(j) < length(block_ids) && block_ids(bounds(j)) ~= block_ids(t - 1)
            break
        end
        theory_id_seq(t:bounds(j)-1) = theory_id_seq(t - 1);

        tic
        nlz_cand = get_theory_seq_nlz(unique_theory_HRRs, theory_id_seq_cand, ts, run_id_frames, R, K, W, x, y, meanfun, covfun, likfun, SPM);
        toc

        fprintf('%d -> %d: nlz = %.5f vs. %.5f\n', t, bounds(j), nlz_cand, nlz_right);
        if nlz_cand <= nlz_right
            fprintf('                        new best: %d vs. %d\n', nlz_cand, nlz_right);
            nlz_right = nlz_cand;
            j_right = j;
            theory_id_seq_right = theory_id_seq_cand;
        else
            disp('                                          not better');
            break
        end

        plot(theory_id_seq);
        hold on;
        plot(block_ids * 100 + 100);
        plot(bounds(j), theory_id_seq(t - 1), '*', 'color', 'red');
        hold off;
        drawnow
        %pause(1);
    end


    % figure out if left or right was better
    if nlz_left < nlz_right
        nlz = nlz_left;
        theory_id_seq = theory_id_seq_left;
        j = j_left;
    else
        nlz = nlz_right;
        theory_id_seq = theory_id_seq_right;
        j = j_right;
    end

    % assess fit
    assert(nlz <= nlz_best); % we're doing greedy => we gotta get better
    if nlz < nlz_best
        fprintf('             new best: %d vs. %d\n', nlz, nlz_best);
        nlz_best = nlz;
        theory_id_seq_best = theory_id_seq;
    end

    % set theory change frame for next interation
    t = bounds(j);
end


theory_id_seq_best_step1 = theory_id_seq_best;
nlz_best_step1 = nlz_best;

save(filename, '-v7.3');



% backwards
%>> scratch
%original
%
%nlz =
%
%   1.8703e+03
%
%
%r =
%
%    0.0408
%
%>> scratch
%best
%
%nlz =
%
%   1.8654e+03
%
%
%r =
%
%    0.1581
%
%>> 


%
% step 2) then fit the actual theories
% use MCMC within Gibbs -- sample from marginal likelihood
% do it stochastically this time b/c now you are changing entire intervals, which might substantially
% change fit for other theories too
%



% for MCMC
theory_id_seq = theory_id_seq_best;
nlz = nlz_best;


tc = find(theory_id_seq(2:end) ~= theory_id_seq(1:end-1)) + 1; % frame indices of theory changes 
tc = [1 tc length(theory_id_seq)+1];

old_bounds = bounds;
bounds = tc; % only consider theory change points
%bounds = cleanup_bounds(bounds); DO NOT cleanup, even though a bit underpowered for some intervals

%while (true) % repeat until we keep improving
for iter = 1:10000
    tic

    % pick random interval
    b = randi(length(bounds)-1);
    %for b = 1:length(bounds) - 1 % for each time point / set of consecutive time points
    st = bounds(b);
    en = bounds(b+1)-1;

    fprintf('between %d and %d\n', st, en);

    % pick theory id candidate at random
    %assert(game_ids(st) == game_ids(en));
    %assert(block_ids(st) == block_ids(en));
    g = game_ids(st);
    tid = randsample(tid_candidates{g}, 1);
    %for tid = 0:n_theories-1 % try assigning each candidate theory to that time point; note b/c of python, tid's are 0-indexed

    fprintf('     assign %d\n', tid);

    % candidate theory sequence = current theory, with one interval changed
    theory_id_seq_cand = theory_id_seq;
    theory_id_seq_cand(st:en) = tid;

    nlz_cand = get_theory_seq_nlz(unique_theory_HRRs, theory_id_seq, ts, run_id_frames, R, K, W, x, y, meanfun, covfun, likfun, SPM);

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
        iter
        save(filename, '-v7.3');
    end

    toc
end

        
save(filename, '-v7.3');
disp('Done');


% temporary (nsamples = 1) to get the other stuff that we didn't log before, but we log now TODO re-run HRR.py w/ nsamples = 100
load(sprintf('mat/unique_HRR_subject_subj=%d_K=10_N=10_E=0.050_nsamples=1_norm=1.mat', subj_id), 'play_key', 'gameStrings', 'unique_theories_filename');

play_key_seq = cell(size(play_key,1),1);
for i = 1:size(play_key_seq,1)
    play_key_seq{i} = play_key(i,:);
end

save(filename, '-v7.3');
disp('Done 2');

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



function nlz = get_theory_seq_nlz(unique_theory_HRRs, theory_id_seq, ts, run_id_frames, R, K, W, x, y, meanfun, covfun, likfun, SPM)

    % get kernel for candidate theory sequence
    [theory_kernel, ~, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq, ts, run_id_frames, SPM);

    ker = R*K*W*theory_kernel*W'*K'*R';

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
    %{
    hyp = hyp_from_ker(nearestSPD(ker));
    hyp.lik = log(1); % TODO sigma_n = 1 const
    nlz = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);
    %}
    nlz = fit_gp_simple(y, ker, x, y, meanfun, covfun, likfun);
    nlz = mean(nlz);

end


function bounds = cleanup_bounds(bounds)

    % optionally remove 1-frame short intervals by concatenating them
    %
    %{
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

    % keep merging smallest interval into neighbor,
    % until there are no small intervals
    % o/w, the temporal resolution is insufficient
    %
    bounds_orig = bounds;
    while true
        intervals = diff(bounds);

        % find smallest interval
        [min_interval, i] = min(intervals);
        if min_interval >= 300 %60 % 1 TR ~= 30 frames
            break
        end

        % get neighbors
        left_interval = inf;
        if i - 1 >= 1
            left_interval = intervals(i - 1);
        end
        right_interval = inf;
        if i + 1 <= length(intervals)
            right_interval = intervals(i + 1);
        end

        % merge w/ shorter neighbor
        if left_interval < right_interval
            bounds(i) = [];
        else
            bounds(i+1) = [];
        end
    end
end
