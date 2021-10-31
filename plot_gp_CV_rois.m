clear all;
close all;

sem = @(x) std(x) / sqrt(length(x));

fasse_ncf = false;
filename = fullfile(get_mat_dir(fasse_ncf), 'gp_CV_rois.mat');
filename

load(filename);

figure;

mean_fs = mean(fs, 3);
h = bar(mean_fs);

xs = [];
for i = 1:nregressors
    xs = [xs, h(i).XData' + h(i).XOffset'];
end

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

ylim([0 0.02])
title('Fraction significant voxels');




figure;

mean_rs = mean(rs, 3);
h = bar(mean_rs);

hold on;
for m = 1:nROIs
    for reg = 1:nregressors
        y = squeeze(rs(m,reg,:));
        x = repmat(xs(m,reg), 1, nsubjects);
        %swarmchart(x, y, 10, h(reg).FaceColor, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5, 'XJitterWidth', 0.05);
        errorbar(xs(m,reg), mean(y), sem(y), '.', 'MarkerSize', 1, 'MarkerFaceColor', h(reg).FaceColor, 'LineWidth', 1, 'Color', h(reg).FaceColor, 'AlignVertexCenters', 'off');
        %Violin(squeeze(fs(m,reg,:)), xs(m,reg), 'Width', 0.1, 'ViolinColor', h(reg).FaceColor, 'ViolinAlpha', 0.3);
    end
end

%ylim([0 0.02])
title('Pearson correlation');
