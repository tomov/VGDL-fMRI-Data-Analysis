clear all;
close all;

sem = @(x) std(x) / sqrt(length(x));

fasse_ncf = false;
filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois.mat');
filename

load(filename);

%% fraction significant voxel z
%

figure;

mean_fs = mean(fs, 3);
h = bar(mean_fs);

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
        %swarmchart(x, y, 10, h(reg).FaceColor, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5, 'XJitterWidth', 0.05);
        %errorbar(xs(m,reg), mean(y), sem(y), '.', 'MarkerSize', 1, 'MarkerFaceColor', h(reg).FaceColor, 'LineWidth', 1, 'Color', h(reg).FaceColor, 'AlignVertexCenters', 'off');
        %Violin(squeeze(fs(m,reg,:)), xs(m,reg), 'Width', 0.1, 'ViolinColor', h(reg).FaceColor, 'ViolinAlpha', 0.3);
    end
end

% significance labels
for m = 1:nROIs
    maxy = max(mean_fs(m,:));
    for r1 = 1:nregressors
        for r2 = r1+1:nregressors
            y1 = squeeze(fs(m,r1,:));
            y2 = squeeze(fs(m,r2,:));
            x1 = xs(m,r1);
            x2 = xs(m,r2);
            %[h,p,ci,stats] = ttest(y1,y2);
            p = ranksum(y1, y2);
            if p <= 0.05
                plot([x1 x2], [maxy maxy] + 0.001, '-', 'color', [0 0 0]);
                text(mean([x1 x2]), maxy + 0.002, significance(p), 'HorizontalAlignment', 'center');
                maxy = maxy + 0.003;
            end

        end
    end
end

legend(regressor_names);
xticklabels(roi_names);
set(gca,'TickLabelInterpreter','none');

%ylim([0 0.06])
title('Fraction significant voxels');


%% Pearson correlations
%

figure;

rs = atanh(rs);

mean_rs = mean(rs, 3);
h = bar(mean_rs);

hold on;
for m = 1:nROIs
    for reg = 1:nregressors
        y = squeeze(rs(m,reg,:));
        x = repmat(xs(m,reg), 1, nsubjects);
        %swarmchart(x, y, 10, h(reg).FaceColor, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5, 'XJitterWidth', 0.05);
        %errorbar(xs(m,reg), mean(y), sem(y), '.', 'MarkerSize', 1, 'MarkerFaceColor', h(reg).FaceColor, 'LineWidth', 1, 'Color', h(reg).FaceColor, 'AlignVertexCenters', 'off');
        %Violin(squeeze(fs(m,reg,:)), xs(m,reg), 'Width', 0.1, 'ViolinColor', h(reg).FaceColor, 'ViolinAlpha', 0.3);
    end
end

% significance labels
for m = 1:nROIs
    maxy = max(mean_rs(m,:));
    for r1 = 1:nregressors
        for r2 = r1+1:nregressors
            y1 = squeeze(rs(m,r1,:));
            y2 = squeeze(rs(m,r2,:));
            x1 = xs(m,r1);
            x2 = xs(m,r2);
            [h,p,ci,stats] = ttest(y1,y2);
            %p = ranksum(y1, y2);
            if p <= 0.05
                plot([x1 x2], [maxy maxy] + 0.001, '-', 'color', [0 0 0]);
                text(mean([x1 x2]), maxy + 0.002, significance(p), 'HorizontalAlignment', 'center');
                maxy = maxy + 0.003;
            end

        end
    end
end

legend(regressor_names);
xticklabels(roi_names);
set(gca,'TickLabelInterpreter','none');

%ylim([0 0.02])
title('Pearson correlation');
