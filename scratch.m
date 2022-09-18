
%clear all;
%close all


% Theory updates
%
%{
EXPT = vgdl_expt;
glmodel = 102;
subj_id = 1;
[~, theory_update] = load_GLM_kernel(EXPT, glmodel, subj_id, {'theory_change_flag'}, false, false);
%}


% Theory representations
%
%{
load(sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', subj_id), 'theory_Xx');
rng default; % reproducibility
Y = tsne(theory_Xx);
%}

[coeff,score,latent,tsquared,explained,mu] = pca(theory_Xx);

figure;
plot(theory_update);
hold on;
plot(Y(:,1));
plot(Y(:,2));
legend({'theory update', 'Y1', 'Y2'});

figure;
scatter(theory_update, Y(:,1));

%for subj = 1:32
%    get_game_for_each_TR(subj, true);
%end


%{
clear all;
close all


EXPT = vgdl_expt;
subj = 1;
[ker, features] = load_state_kernel(EXPT, subj, 1); 

glmodel = 1
[R, K, W] = get_R_K_W(EXPT, glmodel, subj);

figure;
imagesc(ker);
colorbar
figure;
imagesc(features);
colorbar
%}





%{
figure_scale = 0.7;

EXPT = vgdl_expt;
glmodel = 1;
subj_id = 1;

[games, levels] = get_game_for_each_TR(subj_id);

% Get game for each TR
%

% Convert game names

proper_games = convert_game_names(games);

% Get HRR

load(sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', subj_id), 'theory_Xx');

assert(length(games) == size(theory_Xx, 1));
assert(length(levels) == size(theory_Xx, 1));

% Run t-SNE

disp('running tsne...');
tic
rng default; % reproducibility
Y = tsne(theory_Xx);
toc

% Plot t-SNE for all TRs

figure('pos', [712 152 figure_scale*764 figure_scale*764]);
gscatter(Y(:,1), Y(:,2), proper_games);
xlabel('dimension 1');
ylabel('dimension 2');
title('EMPA theory HRRs: t-SNE');
legend('Location','southwest');


% Plot t-SNE for each games

figure('pos', [712 152 figure_scale*764*3/2+300 figure_scale*764]);

game_names_ordered = get_game_names_ordered(subj_id);

proper_game_names_ordered = convert_game_names(game_names_ordered);


for g = 1:6
    subplot(2, 3, g);
    which_rows = strcmp(games, game_names_ordered{g});

    x = Y(which_rows, 1);
    y = Y(which_rows, 2);
    l = (1:sum(which_rows)) / sum(which_rows) * 9;
    hp = patch([x' NaN], [y' NaN], 0);
    set(hp,'cdata', [l NaN], 'edgecolor','interp','facecolor','none');
    hold on;
    scatter(x,y,15,l,'filled');
    hold off;
    hc = colorbar;
    xlabel('dimension 1');
    ylabel('dimension 2');
    ylabel(hc, 'level', 'FontSize', 11);
    title(proper_game_names_ordered{g});
end

%
%disp('running tsne...');
%tic
%rng default; % reproducibility
%[R, K, W] = get_R_K_W(EXPT, glmodel, subj_id);
%Y = tsne(R * K * W * theory_Xx);
%toc
%
%figure;
%gscatter(Y(:,1), Y(:,2), proper_games);
%xlabel('dimension 1');
%ylabel('dimension 2');
%title(sprintf('Subject %d', subj_id));

close(conn);

%}
