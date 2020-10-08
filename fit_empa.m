
clear all;

%load('mat/fmri_empaLik_reg_1_vgfmri3_helper.mat');
%load('mat/fmri_empaLik_best_1_vgfmri3_helper.mat');
load('mat/fmri_empaLik_orig_1_vgfmri3_helper.mat');

data(1).behavior = behavior;
data(1).predictions = predictions;


x = [0.1];

loglik = lik_empa(x, data(1));

loglik
