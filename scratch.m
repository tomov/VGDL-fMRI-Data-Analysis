
figure_scale = 0.7;


%{
% Get game for each TR
%
mongo_connect;

nruns = 6;
nTRs = 1698;
initial_TRs = 7;
n_games_per_run = 3;
nTRs_per_game = (nTRs / nruns - initial_TRs) / n_games_per_run;

EXPT = vgdl_expt;
glmodel = 1;
subj_id = 1;
subj_games = {};

% Get game order for subject
for run_id = 1:nruns
    run = get_run(subj_id, run_id);

    [game_names, onsets, durs] = get_games(subj_id, run, conn);
    subj_games = [subj_games, game_names'];
end

% Expand games and levels to each TR

games = {};
levels = [];

for r = 1:nruns
    games = [games; repmat({''}, [initial_TRs, 1])];
    levels = [levels; repmat(nan, [initial_TRs, 1])];
    
    p = partition_id_from_run_id(r);
    for g = (r - 1) * n_games_per_run + 1 : r * n_games_per_run
        games = [games; repmat(subj_games(g), [nTRs_per_game, 1])];
        levels = [levels; (ceil((1:nTRs_per_game) / nTRs_per_game * 3) + (p - 1) * 3)'];
    end
end

assert(length(games) == nTRs);
assert(length(levels) == nTRs);

% Convert game names

proper_games = games;
proper_games(strcmp(games, 'vgfmri3_chase')) = {'Chase'};
proper_games(strcmp(games, 'vgfmri3_helper')) = {'Helper'};
proper_games(strcmp(games, 'vgfmri3_bait')) = {'Bait'};
proper_games(strcmp(games, 'vgfmri3_lemmings')) = {'Lemmings'};
proper_games(strcmp(games, 'vgfmri3_plaqueAttack')) = {'Plaque Attack'};
proper_games(strcmp(games, 'vgfmri3_zelda')) = {'Zelda'};

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

%}

% Plot t-SNE for each games

figure('pos', [712 152 figure_scale*764*3/2+300 figure_scale*764]);

game_names_ordered = get_game_names_ordered(subj_id);

proper_game_names_ordered = game_names_ordered;
proper_game_names_ordered(strcmp(game_names_ordered, 'vgfmri3_chase')) = {'Chase'};
proper_game_names_ordered(strcmp(game_names_ordered, 'vgfmri3_helper')) = {'Helper'};
proper_game_names_ordered(strcmp(game_names_ordered, 'vgfmri3_bait')) = {'Bait'};
proper_game_names_ordered(strcmp(game_names_ordered, 'vgfmri3_lemmings')) = {'Lemmings'};
proper_game_names_ordered(strcmp(game_names_ordered, 'vgfmri3_plaqueAttack')) = {'Plaque Attack'};
proper_game_names_ordered(strcmp(game_names_ordered, 'vgfmri3_zelda')) = {'Zelda'};


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

%close(conn);


