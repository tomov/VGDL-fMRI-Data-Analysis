function h = plot_gp_CV_rois_helper(fs, test_type, statistic, regressor_names, roi_names, null_value, cmap, significant_scale, regressors_to_compare, font_scale)
    sem = @(x) std(x) / sqrt(length(x));

    nROIs = size(fs, 1);
    nregressors = size(fs, 2);
    nsubjects = size(fs, 3);

    if ~exist('significant_scale', 'var') || isempty(significant_scale)
        significant_scale = 1;
    end
    if ~exist('font_scale', 'var') || isempty(font_scale)
        font_scale = 1;
    end
    if ~exist('regressors_to_compare', 'var') %|| isempty(regressors_to_compare)
        regressors_to_compare = 1:nregressors;
    end

    % bars
    switch statistic
        case 'mean'
            m_fs = mean(fs, 3);
        case 'median'
            m_fs = median(fs, 3);
    end
    sem_fs = std(fs, 0, 3) / sqrt(nsubjects);
    h = bar(1:nROIs, m_fs);

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
                    l_fs(m,reg) = mean(y) - sem(y); % lower confidence bound
                case 'median'
                    q = quantile(y, [0.25 0.75]);
                    neg = m_fs(m,reg) - q(1);
                    pos = q(2) - m_fs(m,reg);
                    errorbar(xs(m,reg), m_fs(m,reg), neg, pos, '.', 'MarkerSize', 1, 'MarkerFaceColor', h(reg).FaceColor, 'LineWidth', 1, 'Color', h(reg).FaceColor, 'AlignVertexCenters', 'off');
                    u_fs(m,reg) = q(2); % upper confidence bound
                    l_fs(m,reg) = q(1); % lower confidence bound
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
                    %text(x1, u_fs(m,r1) + 0.002 * significant_scale, significance(p), 'HorizontalAlignment', 'center');
                    if mean(y1) < 0
                        y = l_fs(m,r1) - 0.003 * significant_scale;
                    else
                        y = u_fs(m,r1) + 0.003 * significant_scale;
                    end
                    text(x1, y, significance(p), 'HorizontalAlignment', 'center', 'FontSize', significant_scale * font_scale);
                   % text(x1, y,  [roi_names{m}, '-', regressor_names{r1}], 'HorizontalAlignment', 'center', 'FontSize', significant_scale);
            %        maxy = maxy + 0.003;
                end
            end
        end

        % paired
        %maxy = max(m_fs(m,:) + sem_fs(m,:));
        maxy = max(u_fs(m,:));
        if exist('null_value') && ~isempty(null_value)
            maxy = maxy + 0.003 * significant_scale; % account for single tests
        end
        for r1 = regressors_to_compare
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
                        %[p,~,stats] = signrank(y1, y2, 'tail', 'right');
                    case 'ranksum'
                        [p,~,stats] = ranksum(y1, y2);
                end
                if p <= 0.05
                    plot([x1 x2], [maxy maxy] + 0.002 * significant_scale, '-', 'color', [0 0 0]);
                    text(mean([x1 x2]), maxy + 0.003 * significant_scale, significance(p), 'HorizontalAlignment', 'center', 'FontSize', significant_scale * font_scale);
                    %text(mean([x1 x2]), maxy + 0.002 * significant_scale, [roi_names{m}, '-', regressor_names{r1}], 'HorizontalAlignment', 'center', 'FontSize', significant_scale);
                    maxy = maxy + 0.003 * significant_scale;
                end
                fprintf('ROI %s: %s vs. %s -- p = %.5f\n', roi_names{m}, regressor_names{r1}, regressor_names{r2}, p);

            end
        end
    end

    legend(regressor_names, 'interpreter', 'none');
    xticks(1:nROIs);
    xticklabels(roi_names);
    xtickangle(30);
    set(gca,'TickLength',[0 0]);
    set(gca,'TickLabelInterpreter','none');

end
