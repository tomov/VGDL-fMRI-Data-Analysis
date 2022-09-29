% aggregate results from fit_gp_CV.m
% copy of agg_gp.m
close all;
clear all;


use_smooth = true;

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

% skip 15 b/c no 3rd partition
subjects = [1:14 16:length(EXPT.subject)];

%subjects = 1:2:32; % odd

%model_name = 'PCA';
%model_name = 'state';
%model_name = 'irrelevant';
%model_name = 'DQN';
%model_name = 'DQN25M';
%model_name = 'game';
model_name = 'EMPA';
%model_name = 'VAE';
%model_name = 'VAE_e1k';
%what = 'conv3';
%%%what = 'linear2';
%what = 'all';
%what = '';
%what = 'novelty';
%what = 'all';
what = 'theory';
%what = 'termination';
%what = 'sprite';
project = 1;
glmodel = 1;
%suffix = '_nowhiten_nofilter';
%suffix = '_';
suffix = '_parts=123';
normalize = 1;
concat = 0;
novelty = 1;
saveYhat = 0;

agg_filename = fullfile(get_mat_dir(), sprintf('agg_gp_CV_us=%d_glm=%d_model=%s_%s_nsamples=100_project=%d_norm=%d_concat=%d_novelty=%d_fast=1%s_by_game.mat', use_smooth, glmodel, model_name, what, project, normalize, concat, novelty, suffix));
agg_filename

for s = 1:1:length(subjects)
    subj_id = subjects(s);

    game_names = get_game_names_ordered(subj_id);
    for g = 1:length(game_names)
        game_name = game_names{g};

        filename = sprintf('/n/holyscratch01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/fit_gp_CV_subj=%d_us=1_glm=1_mask=mask_model=EMPA_theory_nsamples=100_project=1_norm=1_concat=0_novelty=1_fast=1_saveYhat=1_parts=123_games=%s.mat', subj_id, game_name);
        filename

        load(filename, 'n', 'logmarglik', 'adjR2', 'R2_CV', 'sigma', 'mask', 'r', 'r_CV', 'mask', 'MSE', 'SMSE', 'MSE_CV', 'SMSE_CV', 'partition_id');

        if s == 1
            rs = nan(length(subjects), length(game_names), size(r,2));
        end

        rs(s,g,:) = mean(r_CV, 1);
    end

end

% Fisher z transform
zs = atanh(rs);

zs_collapsed = sum(rs, 2);

% t-test Pearson corr across subjects
[h,p,ci,stats] = ttest(zs_collapsed);
ts = stats.tstat;


agg_filename
save(agg_filename);

