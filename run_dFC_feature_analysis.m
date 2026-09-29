function run_dFC_feature_analysis(cfg)
% RUN_DFC_FEATURE_ANALYSIS: Extract features, perform ANOVA/FDR tests and generate graphs.
% Usage：cfg = config_analysis(); run_dFC_feature_analysis(cfg);

    if nargin < 1
        cfg = config_analysis();
    end
    rng(cfg.randomSeed);   
    if ~exist(cfg.resultsRoot, 'dir'), mkdir(cfg.resultsRoot); end
    if ~exist(cfg.figureDir, 'dir'),   mkdir(cfg.figureDir);   end

    nBands    = size(cfg.freqBands, 1);
    nFeatures = 4;

    feature_control = cell(1, nBands);
    feature_good    = cell(1, nBands);
    feature_bad     = cell(1, nBands);

    groups = {cfg.groups.control, cfg.groups.good, cfg.groups.bad};

    % 
    for g = 1:numel(groups)
        for fre = 1:nBands
            data = [];
            for p = 1:numel(groups{g})
                patient = groups{g}{p};
                f = fullfile(cfg.derivRoot, patient, 'wpli_dynamic.mat');
                if ~isfile(f)
                    warning('Cannot find %s, skipping.', f);
                    continue;
                end
                S = load(f, 'wpli_dynamic');
                features = extract_dFC_features(S.wpli_dynamic{1,fre}, cfg.nClusters);
                data = [data; features]; 
            end
            switch g
                case 1, feature_control{fre} = data;
                case 2, feature_good{fre}    = data;
                case 3, feature_bad{fre}     = data;
            end
        end
    end

    % 
    significant_FDR_value = cell(nBands, nFeatures);
    significant_value     = cell(nBands, nFeatures);
    meanmedian_value      = cell(nBands, nFeatures);
    effect_value          = cell(nBands, nFeatures);
    effect_type           = cell(nBands, nFeatures);

    for fre = 1:nBands
        for fea = 1:nFeatures
            data1 = feature_control{fre}(:, fea);
            data2 = feature_good{fre}(:, fea);
            data3 = feature_bad{fre}(:, fea);

            [significant_value1, meanmedian_value1, effect_table1] = ...
                ANOVA_effect_size(data1, data2, data3);
            close all;

            if size(significant_value1, 2) > 3
                significant_value2 = abs(significant_value1(:,6))';
            else
                significant_value2 = [significant_value1(1,2), ...
                    significant_value1(1,3), significant_value1(2,3)];
            end
            fdr_adj = mafdr(significant_value2, 'BHFDR', true);

            significant_FDR_value{fre,fea} = fdr_adj;
            significant_value{fre,fea}     = significant_value2;
            meanmedian_value{fre,fea}      = meanmedian_value1;
            effect_value{fre,fea}          = (effect_table1.EffectSize)';

            type1 = [];
            for ii = 1:3
                if strcmp(effect_table1.Type{ii}, "Hedges' g")
                    type11 = 1;
                elseif strcmp(effect_table1.Type{ii}, "Glass's delta (vs HC SD)")
                    type11 = 2;
                else
                    type11 = 3;
                end
                type1 = [type1, type11]; 
            end
            effect_type{fre,fea} = type1;
        end
    end

    % 
    save(fullfile(cfg.resultsRoot, 'dFC_statistics.mat'), ...
        'significant_FDR_value', 'significant_value', 'meanmedian_value', ...
        'effect_value', 'effect_type', 'feature_control', 'feature_good', ...
        'feature_bad', '-v7.3');

    % 
    fig = figure('Units','centimeters','Position',[0.1 0.1 48 25], ...
        'Visible','off');
    colors = [100 149 237; 254 194 140; 124 208 159] / 255;
    jitterWidth = 0.2;

    for fre = 1:nBands
        data1 = feature_control{fre}(:,1);
        data2 = feature_good{fre}(:,1);
        data3 = feature_bad{fre}(:,1);

        means = [mean(data1), mean(data2), mean(data3)];
        sems  = [std(data1), std(data2), std(data3)];
        dataCell = {data1, data2, data3};
        nGroups = numel(dataCell);

        subplot(2,3,fre);
        b = bar(1:3, means, 'FaceColor','flat', ...
            'EdgeColor','flat', 'BarWidth',0.56);
        b.CData = colors;
        b.FaceAlpha = 0.3;
        b.LineWidth = 3;
        hold on;

        set(gca, 'XTickLabel', cfg.groupLabels, ...
            'FontSize', 30, 'FontWeight','bold');
        ylabel('Global Variability', 'FontSize', 33, 'FontWeight','bold');
        set(gca, 'LineWidth', 5, 'FontWeight','bold', 'FontSize', 25);
        box off;

        try
            xCenters = b.XEndPoints;
        catch
            xCenters = 1:nGroups;
        end

        for i = 1:nGroups
            y = dataCell{i};
            xPos = xCenters(i);
            xJitter = xPos + (rand(size(y)) - 0.5) * jitterWidth;
            scatter(xJitter, y, 40, 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', colors(i,:));
        end

        for i = 1:nGroups
            errorbar(xCenters(i), means(i), NaN, sems(i), ...
                'k', 'LineWidth', 3, 'CapSize', 17);
        end
    end

    exportgraphics(fig, fullfile(cfg.figureDir, ...
        'dFC_GlobalVariability.png'), 'Resolution', 300);
    close(fig);
end
