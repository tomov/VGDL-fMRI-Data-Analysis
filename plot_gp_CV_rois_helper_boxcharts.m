function h = plot_gp_CV_rois_helper_boxcharts(fs, test_type, statistic, regressor_names, roi_names, null_value, cmap, significant_scale, regressors_to_compare, font_scale, outliers)
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
    if ~exist('outliers', 'var') || isempty(outliers)
        outliers = false;
    end
    xscale = 1.25;

    % convert to table, for box charts
    x_rows = [];
    roi_rows = [];
    regressor_rows = [];
    subject_rows = [];
    f_rows = [];
    for m = 1:nROIs
        for reg = 1:nregressors
            y = squeeze(fs(m,reg,:));
            %x = m - 0.2 + reg * 0.1;
            x = m * xscale;
            x_rows = [x_rows; repmat(x, length(y), 1)];
            roi_rows = [roi_rows; repmat(roi_names(m), length(y), 1)];
            regressor_rows = [regressor_rows; repmat(regressor_names(reg), length(y), 1)];
            subject_rows = [subject_rows; (1:length(y))'];
            f_rows = [f_rows; y];
            if outliers
                u_fs(m,reg) = max(y); % upper limit
            else
                u_fs(m,reg) = max(y(~isoutlier(y, 'quartiles'))); % upper limit
            end
            xs(m,reg) = x - 0.45 / 1.5 * xscale + 0.75 * (reg - 1) / (nregressors - 1);
        end
    end
    tbl = table(x_rows, roi_rows, regressor_rows, subject_rows, f_rows, 'VariableNames', {'x', 'ROI', 'regressor', 'subject', 'f'});
    tbl.ROI = categorical(tbl.ROI, roi_names);
    tbl.regressor = categorical(tbl.regressor, regressor_names);

    % bars
    switch statistic
        case 'mean'
            assert(false, 'Box charts only make sense for medians');
        case 'median'
            m_fs = median(fs, 3);
            %b = boxchart(tbl.ROI, tbl.f, 'GroupByColor', tbl.regressor, 'Notch', 'on', 'BoxWidth', 0.75);
            h = boxchart(tbl.x, tbl.f, 'GroupByColor', tbl.regressor, 'BoxWidth', 0.75); 
            if outliers
                for i = 1:nregressors
                    h(i).JitterOutliers = 'on';
                    h(i).MarkerStyle = '.';
                end
            else
                for i = 1:nregressors
                    h(i).MarkerStyle = 'none';
                end
            end
    end

    % color the bars
    if exist('cmap', 'var') && ~isempty(cmap)
        for i = 1:length(h)
            h(i).BoxFaceColor = cmap(i,:);
            h(i).MarkerColor = cmap(i,:);
        end
    end


    hold on;

    % significance labels
    for m = 1:nROIs

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
                    case 'ranksum'
                        [p,~,stats] = ranksum(y1, y2);
                end
                if p <= 0.05
                    plot([x1 x2], [maxy maxy] + 0.002 * significant_scale, '-', 'color', [0 0 0]);
                    text(mean([x1 x2]), maxy + 0.0025 * significant_scale, significance(p), 'HorizontalAlignment', 'center', 'FontSize', significant_scale * font_scale);
                    %text(mean([x1 x2]), maxy + 0.002 * significant_scale, [roi_names{m}, '-', regressor_names{r1}], 'HorizontalAlignment', 'center', 'FontSize', significant_scale);
                    maxy = maxy + 0.003 * significant_scale;
                end
                fprintf('ROI %s: %s vs. %s -- p = %.5f\n', roi_names{m}, regressor_names{r1}, regressor_names{r2}, p);

            end
        end
    end

    legend(regressor_names, 'interpreter', 'none');
    xticks(xscale:xscale:xscale*nROIs);
    xticklabels(roi_names);
    xtickangle(30);
    set(gca,'TickLength',[0 0]);
    set(gca,'TickLabelInterpreter','none');

end
