clear all;
close all;


fasse_ncf = false;
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=HarvardOxford-maxprob-thr0.mat');
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL3v1.mat');

agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2.mat'); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010_atlas=AAL2_grouped.mat');

%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.010.mat');
%agg_filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois_alpha=0.001.mat');
agg_filename

load(agg_filename);

%
%% BICs
%

%{
figure('position', [673 90 1519 849]);
ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
h = plot_gp_CV_rois_helper(bics(:,ix,:) - bics(:,end,:), 'ttest', 'mean', regressor_names(ix), roi_names);
title('BICs in ROIs');
ylabel('\Delta BIC');
%}

%
%% fraction significant voxel z
%

figure('position', [673 90 1519 849]);
ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names);
%h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, alpha);
%plot_gp_CV_rois_helper(fs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names);
%ylim([0 0.1]);
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');


%{
figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'theory', 'sprite', 'interaction', 'termination'});
cmap = [1 0.8 0.6 0.4]' * h(1).FaceColor;
plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, alpha); %colormap(winter(3)));
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');

figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'DQN', 'conv1', 'conv2', 'conv3', 'linear1', 'linear2'});
cmap = [1 0.9 0.8 0.7 0.6 0.5]' * h(2).FaceColor;
h = plot_gp_CV_rois_helper(fs(:,ix,:), 'signrank', 'median', regressor_names(ix), roi_names, alpha); %colormap(autumn(5)));
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');
%}


%{
figure;
plot_gp_CV_rois_helper(fs, 'signrank', regressor_names, roi_names);
title('Fraction significant voxels in ROIs');
ylabel('Fraction significant voxels');
%}

%% Pearson correlations
%


%{
figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'theory', 'DQN', 'PCA'});
plot_gp_CV_rois_helper(zs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names, 0);
title('Fisher z-transformed Pearson correlation between predicted and actual BOLD');
ylabel('Fisher z-transformed Pearson correlation coefficient');
%}


%{
figure;
plot_gp_CV_rois_helper(zs, 'ttest', regressor_names, roi_names);
title('Pearson correlation');
%}

%{
figure('position', [73 90 1519 849]);
ix = ismember(regressor_names, {'sprite', 'interaction', 'termination'});
h = plot_gp_CV_rois_helper(zs(:,ix,:), 'ttest', 'mean', regressor_names(ix), roi_names);
%ylim([0 0.1]);
title('Fisher z-transformed Pearson correlation between predicted and actual BOLD');
ylabel('Fisher z-transformed Pearson correlation coefficient');
cmap = colormap(winter(3));
for i = 1:length(h)
    h(i).FaceColor = cmap(i,:);
end
%}



function h = plot_gp_CV_rois_helper(fs, test_type, statistic, regressor_names, roi_names, null_value, cmap)
    sem = @(x) std(x) / sqrt(length(x));

    nROIs = size(fs, 1);
    nregressors = size(fs, 2);
    nsubjects = size(fs, 3);

    % bars
    switch statistic
        case 'mean'
            m_fs = mean(fs, 3);
        case 'median'
            m_fs = median(fs, 3);
    end
    sem_fs = std(fs, 0, 3) / sqrt(nsubjects);
    h = bar(m_fs);

    % color the bars
    if exist('cmap', 'var') && ~isempty(cmap)
        for i = 1:length(h)
            h(i).FaceColor = cmap(i,:);
        end
    end

    xs = [];
    for i = 1:nregressors
        xs = [xs, h(i).XData' + h(i).XOffset'];
    end

    % error bars
    hold on;
    for m = 1:nROIs
        for reg = 1:nregressors
            y = squeeze(fs(m,reg,:));
            x = repmat(xs(m,reg), 1, nsubjects);
            u_fs(m,reg) = mean(y); % default upper confidence bound
            %swarmchart(x, y, 10, h(reg).FaceColor, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5, 'XJitterWidth', 0.15);
            switch statistic
                case 'mean'
                    errorbar(xs(m,reg), mean(y), sem(y), '.', 'MarkerSize', 1, 'MarkerFaceColor', h(reg).FaceColor, 'LineWidth', 1, 'Color', h(reg).FaceColor, 'AlignVertexCenters', 'off');
                    u_fs(m,reg) = mean(y) + sem(y); % upper confidence bound
                case 'median'
                    q = quantile(y, [0.25 0.75]);
                    neg = m_fs(m,reg) - q(1);
                    pos = q(2) - m_fs(m,reg);
                    errorbar(xs(m,reg), m_fs(m,reg), neg, pos, '.', 'MarkerSize', 1, 'MarkerFaceColor', h(reg).FaceColor, 'LineWidth', 1, 'Color', h(reg).FaceColor, 'AlignVertexCenters', 'off');
                    u_fs(m,reg) = q(2); % upper confidence bound
            end
            %Violin(squeeze(fs(m,reg,:)), xs(m,reg), 'Width', 0.1, 'ViolinColor', h(reg).FaceColor, 'ViolinAlpha', 0.3);
        end
    end

    % significance labels
    for m = 1:nROIs

        % single
        if exist('null_value') && ~isempty(null_value)
            for r1 = 1:nregressors
                y1 = squeeze(fs(m,r1,:));
                x1 = xs(m,r1);
                switch test_type
                    case 'ttest'
                        [~,p,ci,stats] = ttest(y1, null_value);
                    case 'signrank'
                        [p,~,stats] = signrank(y1, null_value, 'Tail','right');
                end
                if p <= 0.05
                    text(x1, u_fs(m,r1) + 0.002, significance(p), 'HorizontalAlignment', 'center');
            %        maxy = maxy + 0.003;
                end
            end
        end

        % paired
        %maxy = max(m_fs(m,:) + sem_fs(m,:));
        maxy = max(u_fs(m,:));
        if exist('null_value') && ~isempty(null_value)
            maxy = maxy + 0.003; % account for single tests
        end
        for r1 = 1:nregressors
            for r2 = r1+1:nregressors
                y1 = squeeze(fs(m,r1,:));
                y2 = squeeze(fs(m,r2,:));
                x1 = xs(m,r1);
                x2 = xs(m,r2);
                switch test_type
                    case 'ttest'
                        [~,p,ci,stats] = ttest(y1,y2);
                    case 'signrank'
                        [p,~,stats] = signrank(y1, y2);
                end
                if p <= 0.05
                    plot([x1 x2], [maxy maxy] + 0.001, '-', 'color', [0 0 0]);
                    text(mean([x1 x2]), maxy + 0.002, significance(p), 'HorizontalAlignment', 'center');
                    maxy = maxy + 0.003;
                end

            end
        end
    end

    legend(regressor_names);
    xticks(1:nROIs);
    xticklabels(roi_names);
    xtickangle(30);
    set(gca,'TickLength',[0 0]);
    set(gca,'TickLabelInterpreter','none');

end
