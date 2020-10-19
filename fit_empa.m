
%clear all;


%load('mat/fmri_empaLik_reg_1_vgfmri3_helper_old.mat');
%load('mat/fmri_empaLik_best_1_vgfmri3_helper_old.mat');
%load('mat/fmri_empaLik_orig_1_vgfmri3_helper_old.mat');

%load('mat/fmri_empaLik_orig_1_vgfmri3_chase.mat');
%load('mat/fmri_empaLik_orig_1_vgfmri3_bait.mat');
%load('mat/fmri_empaLik_orig_1_vgfmri3_helper.mat');
%load('mat/fmri_empaLik_orig_1_vgfmri3_zelda.mat');
%load('mat/fmri_empaLik_orig_1_vgfmri3_lemmings.mat');
%load('mat/fmri_empaLik_orig_1_vgfmri3_plaqueAttack.mat');

%helper(behavior, predictions);


%load('mat/fmri_empaLik_best_1_vgfmri3_chase.mat');
%load('mat/fmri_empaLik_best_1_vgfmri3_bait.mat');
%load('mat/fmri_empaLik_best_1_vgfmri3_helper.mat');
%load('mat/fmri_empaLik_best_1_vgfmri3_zelda.mat');
%load('mat/fmri_empaLik_best_1_vgfmri3_lemmings.mat');
%load('mat/fmri_empaLik_best_1_vgfmri3_plaqueAttack.mat');

%helper(behavior, predictions);

game_names = { ...
    'vgfmri3_chase'; ...
    'vgfmri3_helper'; ...
    'vgfmri3_bait'; ...
    'vgfmri3_zelda'; ...
    'vgfmri3_lemmings'; ...
    'vgfmri3_plaqueAttack'};

original_logliks = [];
decoded_logliks = [];
p_values = [];
LR_stats = [];
original_BICs = [];
decoded_BICs = [];
BIC_diffs = [];

for g = 1:length(game_names)
    [lo, lb] = compare(game_names{g});

    original_logliks = [original_logliks; lo];
    decoded_logliks = [decoded_logliks; lb];
    
    [h,p,stat,c] = lratiotest(lb, lo, 1);
    p_values = [p_values; p];
    LR_stats = [LR_stats; stat];
    
    original_BICs = [original_BICs; - 2 * lo];
    decoded_BICs = [decoded_BICs; - 2 * lb];
    BIC_diffs = [BIC_diffs; -2*lo - (-2*lb)];
end

T = table(game_names, original_logliks, decoded_logliks, p_values, LR_stats, original_BICs, decoded_BICs, BIC_diffs)



% log lik original, log lik best (decoded)
function [lo, lb] = compare(game_name)

    load(sprintf('mat/fmri_empaLik_orig_1_%s.mat', game_name));
    pko = play_keys;
    if ischar(pko)
        pko = cellstr(pko);
    end

    load(sprintf('mat/fmri_empaLik_best_1_%s.mat', game_name));
    pkb = play_keys;
    if ischar(pkb)
        pkb = cellstr(pkb);
    end
   
    l = min(length(pko), length(pkb));
    pko = pko(1:l);
    pkb = pkb(1:l);
    which = strcmp(pko, pkb);

    
    load(sprintf('mat/fmri_empaLik_orig_1_%s.mat', game_name));
    [orig, lo] = helper(behavior, predictions, which);

    load(sprintf('mat/fmri_empaLik_best_1_%s.mat', game_name));
    [best, lb] = helper(behavior, predictions, which);

end



function [data, loglik] = helper(behavior, predictions, which)

    x = [0.1];
    if ischar(behavior)
        data(1).behavior = cellstr(squeeze(behavior));
    else
        data(1).behavior = behavior;
    end
    data(1).predictions = predictions;

    if ~exist('which', 'var')
        which = logical(ones(1, length(data(1).behavior)));
    end
    data(1).behavior = data(1).behavior(which);
    data(1).predictions = data(1).predictions(which);

    loglik = lik_empa(x, data(1));
end
