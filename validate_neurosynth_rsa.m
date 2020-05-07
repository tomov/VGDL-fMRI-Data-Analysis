clear all;

load('mat/neurosynth_rsa_2_us=0_l=1_nperms=100_nroi=7_pi=177,210,352,143,106,38,389.mat');

if use_smooth
    EXPT = vgdl_expt();
else
    EXPT = vgdl_expt_nosmooth();
end

subjects = 1:8;

names = region(which);

[Neural] = ccnl_roi_rdms(EXPT, rsa_idx, roi_masks, subjects);

%% from vgdl_create_rsa
game_names = {'vgfmri3_chase','vgfmri3_helper','vgfmri3_bait','vgfmri3_lemmings','vgfmri3_plaqueAttack','vgfmri3_zelda'};
game_name_to_id = containers.Map(game_names, 1:6);

m = 1; % model

close all;

for i = 1:length(roi_masks)

    parcel_idx = names(i);

    %clear cross_dist, within_dist;

    for s = 1:length(Neural(i).subj) % subject
        upper = logical(triu(ones(size(Behavioral(m).subj(s).RDM)), 1));
        [rx,ry] = meshgrid(Behavioral(m).subj(s).runs);
        [gx,gy] = meshgrid(Behavioral(m).subj(s).features);

        for r = 1:3 % run pair
            cross_game = gx ~= gy & upper & rx == r & ry == r;
            cross_dist(s,r) = mean(Neural(i).subj(s).RDM(cross_game));

            for g = 1:6 % game
                within_game = gx == g & gy == g & upper & rx == r & ry == r;
                within_dist(s,r,g) = mean(Neural(i).subj(s).RDM(within_game));
            end
        end
    end

    % collapse across subjects
    wm = squeeze(mean(within_dist, 1));
    wse = squeeze(std(within_dist, 0, 1)) / sqrt(length(subjects));

    figure;

    % within game, split by game
    %

    subplot(2,3,1);
    errorbar(wm, wse, 'linewidth', 2);
    xlim([0.5 3.5]);
    legend(game_names, 'interpreter', 'none');
    ylabel('avg cosine dist');
    xlabel('run pair');
    title([ num2str(parcel_idx), ': within-game, by game']);

    subplot(2,3,4);
    errorbar(mean(wm',1), std(wm', 0, 1) / sqrt(6), 'linewidth', 2); % TODO proper within-subject SEMs
    xlim([0.5 3.5]);
    ylabel('avg cosine dist');
    xlabel('run pair');
    title([ num2str(parcel_idx), ': within-game']);

    % within game, split by subject
    %

    % collapse across games
    wm_g = squeeze(mean(within_dist, 3));
    wse_g = squeeze(std(within_dist, 0, 3)) / sqrt(6);

    subplot(2,3,2);
    errorbar(wm_g', wse_g', 'linewidth', 2);
    xlim([0.5 3.5]);
    legend(cellfun(@num2str, num2cell(subjects), 'UniformOutput', false), 'interpreter', 'none');
    ylabel('avg cosine dist');
    xlabel('run pair');
    title([ num2str(parcel_idx), ': within-game, by subject']);

    subplot(2,3,5);
    errorbar(mean(wm_g,1), std(wm_g, 0, 1) / sqrt(length(subjects)), 'linewidth', 2); % TODO proper within-subject SEMs
    xlim([0.5 3.5]);
    ylabel('avg cosine dist');
    xlabel('run pair');
    title([ num2str(parcel_idx), ': within-game']);

    % cross game
    %

    subplot(2,3,3);
    plot(cross_dist', 'linewidth', 2);
    xlim([0.5 3.5]);
    legend(cellfun(@num2str, num2cell(subjects), 'UniformOutput', false), 'interpreter', 'none');
    ylabel('avg cosine dist');
    xlabel('run pair');
    title([num2str(parcel_idx), ': cross-game, by subject']);

    subplot(2,3,6);
    errorbar(mean(cross_dist,1), std(cross_dist, 0, 1) / sqrt(numel(subjects)), 'linewidth', 2);
    xlim([0.5 3.5]);
    ylabel('avg cosine dist');
    xlabel('run pair');
    title([num2str(parcel_idx), ': cross-game, by subject']);
end
