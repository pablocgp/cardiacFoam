function renderPETScConvergenceBlocks(dimTag, resultRoot, nRangeAll, nRangeFine)
% Render six-panel convergence blocks for PETSc CN/consistent MMS runs.
%
% Figures are split by nonlinear method, mesh family, stabilisation alpha,
% and time step. Each panel includes inline model labels such as NO-NO and
% p3-p2, following the style of the older JfNK report plots.

outputDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

methodTags = {'Picard', 'diagonalIion', 'JFNK'};
meshTags = {'hexa', 'triangular_Unstr'};
alphaTags = {'0p0e00', '1p0em01'};
dtTags = {'1p0em02', '1p0em03'};
schemeTag = 'CN';
massTag = 'cons';

metrics = struct( ...
    'column', {'Vm_L2', 'VmG_L2', 'u1_L2', 'u1G_L2', 'u2_L2', 'u2G_L2'}, ...
    'label', {'V_m L_2', 'V_m^G L_2', 'u_1 L_2', 'u_1^G L_2', 'u_2 L_2', 'u_2^G L_2'} ...
);

for methodI = 1:numel(methodTags)
    for meshI = 1:numel(meshTags)
        for alphaI = 1:numel(alphaTags)
            for dtI = 1:numel(dtTags)
                curves = collectCurves( ...
                    dimTag, resultRoot, meshTags{meshI}, alphaTags{alphaI}, ...
                    dtTags{dtI}, methodTags{methodI}, schemeTag, massTag, ...
                    metrics, nRangeAll, nRangeFine ...
                );

                if isempty(curves)
                    warning( ...
                        'No curves found for %s method=%s mesh=%s alpha=%s dt=%s', ...
                        dimTag, methodTags{methodI}, meshTags{meshI}, alphaTags{alphaI}, dtTags{dtI} ...
                    );
                    continue;
                end

                renderBlock( ...
                    curves, dimTag, meshTags{meshI}, alphaTags{alphaI}, dtTags{dtI}, ...
                    methodTags{methodI}, schemeTag, massTag, metrics, outputDir ...
                );
            end
        end
    end
end
end

function curves = collectCurves(dimTag, resultRoot, meshTag, alphaTag, dtTag, methodTag, schemeTag, massTag, metrics, nRangeAll, nRangeFine)
curves = struct( ...
    'vm', {}, 'iion', {}, 'N', {}, 'h', {}, 'values', {}, ...
    'ordersAll', {}, 'ordersFine', {}, 'totalTime', {}, 'totalRss', {} ...
);

meshDir = fullfile(resultRoot, sprintf('%s_mesh_%s_alpha_%s', dimTag, meshTag, alphaTag));
if ~exist(meshDir, 'dir')
    return;
end

configPattern = sprintf('hoVm_*_hoStates_*_%s_%s_%s_RKF45', schemeTag, massTag, methodTag);
configDirs = dir(fullfile(meshDir, configPattern));
configRegex = sprintf( ...
    '^hoVm_([^_]+)_hoStates_([^_]+)_%s_%s_%s_RKF45$', ...
    schemeTag, massTag, methodTag ...
);

for configI = 1:numel(configDirs)
    configDir = configDirs(configI);
    if ~configDir.isdir
        continue;
    end

    configMatch = regexp(configDir.name, configRegex, 'tokens', 'once');
    if isempty(configMatch)
        continue;
    end

    endDirs = dir(fullfile(configDir.folder, configDir.name, ['endT_*_dt_' dtTag]));
    if isempty(endDirs)
        continue;
    end

    endDir = endDirs(1);
    errorFile = fullfile(endDir.folder, endDir.name, ['errors_vs_N_dt_' dtTag '.dat']);
    if ~exist(errorFile, 'file')
        continue;
    end

    rows = readErrorFile(errorFile);
    if isempty(rows)
        continue;
    end

    N = meshN(rows);
    [N, order] = sort(N);
    rows = rows(order, :);
    h = 1 ./ N;

    values = nan(height(rows), numel(metrics));
    ordersAll = nan(1, numel(metrics));
    ordersFine = nan(1, numel(metrics));
    for metricI = 1:numel(metrics)
        columnName = metrics(metricI).column;
        if ismember(columnName, rows.Properties.VariableNames)
            values(:, metricI) = rows.(columnName);
            ordersAll(metricI) = fitOrder(N, rows.(columnName), nRangeAll);
            ordersFine(metricI) = fitOrder(N, rows.(columnName), nRangeFine);
        end
    end

    [totalTime, totalRss] = accumulatedCost(errorFile, rows);
    curves(end + 1) = struct( ...
        'vm', configMatch{1}, ...
        'iion', configMatch{2}, ...
        'N', N, ...
        'h', h, ...
        'values', values, ...
        'ordersAll', ordersAll, ...
        'ordersFine', ordersFine, ...
        'totalTime', totalTime, ...
        'totalRss', totalRss ...
    ); %#ok<AGROW>
end

curves = sortCurves(curves);
end

function renderBlock(curves, dimTag, meshTag, alphaTag, dtTag, methodTag, schemeTag, massTag, metrics, outputDir)
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [60 60 1900 1180]);
layout = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for metricI = 1:numel(metrics)
    ax = nexttile(layout);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    legendHandles = gobjects(0);
    legendLabels = {};
    inlineH = {};
    inlineE = {};
    inlineLabels = {};
    inlineColors = zeros(numel(curves), 3);

    for curveI = 1:numel(curves)
        e = curves(curveI).values(:, metricI);
        valid = isfinite(curves(curveI).h) & isfinite(e) & curves(curveI).h > 0 & e > 0;
        if ~any(valid)
            continue;
        end

        curveColor = vmColor(curves(curveI).vm);
        curveMarker = iionMarker(curves(curveI).iion);
        hp = loglog( ...
            ax, curves(curveI).h(valid), e(valid), '-', ...
            'Color', curveColor, ...
            'Marker', curveMarker, ...
            'LineWidth', 1.12, ...
            'MarkerSize', 4.0 ...
        );

        legendHandles(end + 1) = hp; %#ok<AGROW>
        legendLabels{end + 1} = legendCurveLabel(curves(curveI), metricI); %#ok<AGROW>
        inlineH{end + 1} = curves(curveI).h(valid); %#ok<AGROW>
        inlineE{end + 1} = e(valid); %#ok<AGROW>
        inlineLabels{end + 1} = modelLabel(curves(curveI)); %#ok<AGROW>
        inlineColors(numel(inlineLabels), :) = curveColor;
    end

    occupiedBoxes = [];
    for labelI = 1:numel(inlineLabels)
        occupiedBoxes = addInlineCurveLabel( ...
            ax, inlineH{labelI}, inlineE{labelI}, inlineLabels{labelI}, ...
            inlineColors(labelI, :), labelI, occupiedBoxes ...
        );
    end

    xlabel(ax, 'h = 1/N');
    ylabel(ax, metrics(metricI).label, 'Interpreter', 'tex');
    title(ax, metrics(metricI).label, 'Interpreter', 'tex');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');

    if ~isempty(legendHandles)
        lgd = legend(ax, legendHandles, legendLabels, ...
            'Location', 'southwest', 'Interpreter', 'none', 'FontSize', 5.8);
        lgdTitle = title(lgd, 'model | OC all/fine | t_\Sigma | RSS_\Sigma');
        lgdTitle.Interpreter = 'tex';
        if isprop(lgd, 'ItemTokenSize')
            lgd.ItemTokenSize = [9 5];
        end
        lgd.Box = 'on';
        lgd.Color = 'white';
        lgd.EdgeColor = [0.75 0.75 0.75];
    end
end

blockTitle = sprintf( ...
    '%s %s, %s/%s, mesh=%s, alpha=%s, dt=%s', ...
    dimTag, methodTag, schemeTag, massTag, meshLabel(meshTag), alphaLabel(alphaTag), dtLabel(dtTag) ...
);
try
    title(layout, blockTitle, 'Interpreter', 'none');
catch
    sgtitle(blockTitle, 'Interpreter', 'none');
end

baseName = sprintf( ...
    'implicit_%s_%s_%s_%s_%s_alpha_%s_dt_%s', ...
    dimTag, methodTag, schemeTag, massTag, meshTag, alphaTag, dtTag ...
);
pdfFile = fullfile(outputDir, [baseName '.pdf']);
pngFile = fullfile(outputDir, [baseName '.png']);
saveFigure(fig, pdfFile, pngFile);
close(fig);

fprintf('Wrote %s and %s\n', pdfFile, pngFile);
end

function rows = readErrorFile(errorFile)
fid = fopen(errorFile, 'r');
if fid < 0
    rows = table();
    return;
end
cleanupObj = onCleanup(@() fclose(fid));

headerLine = fgetl(fid);
if ~ischar(headerLine)
    rows = table();
    return;
end
rawNames = strsplit(regexprep(strtrim(headerLine), '^#\s*', ''));
varNames = matlab.lang.makeValidName(rawNames);
formatSpec = strjoin(repmat({'%s'}, 1, numel(varNames)));
columns = textscan(fid, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'CommentStyle', '#');

if isempty(columns) || isempty(columns{1})
    rows = table();
    return;
end

rows = table(columns{:}, 'VariableNames', varNames);
textColumns = {'dim', 'mesh_type', 'config', 'dat_file'};
for colI = 1:numel(varNames)
    if ~ismember(varNames{colI}, textColumns)
        rows.(varNames{colI}) = str2double(rows.(varNames{colI}));
    end
end
end

function N = meshN(rows)
N = double(rows.N);
if ismember('dat_file', rows.Properties.VariableNames)
    for rowI = 1:height(rows)
        token = regexp(rows.dat_file{rowI}, '_(\d+)_cells', 'tokens', 'once');
        if ~isempty(token)
            parsedN = str2double(token{1});
            if isfinite(parsedN) && parsedN > 0
                N(rowI) = parsedN;
            end
        end
    end
end
end

function order = fitOrder(N, errorValues, nRange)
valid = isfinite(N) & isfinite(errorValues) & errorValues > 0 & N >= nRange(1) & N <= nRange(2);
if nnz(valid) < 2
    order = NaN;
    return;
end

x = log(1 ./ N(valid));
y = log(errorValues(valid));
xMean = mean(x);
yMean = mean(y);
den = sum((x - xMean).^2);
if den <= 0
    order = NaN;
else
    order = sum((x - xMean) .* (y - yMean)) / den;
end
end

function [totalTime, totalRss] = accumulatedCost(errorFile, rows)
totalTime = 0.0;
totalRss = 0.0;
timeCount = 0;
rssCount = 0;

for rowI = 1:height(rows)
    datFile = fullfile(fileparts(errorFile), rows.dat_file{rowI});
    if ~exist(datFile, 'file')
        continue;
    end

    text = fileread(datFile);
    timeToken = regexp(text, '\<total\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
    rssToken = regexp(text, 'peakRSS_MB\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
    if ~isempty(timeToken)
        totalTime = totalTime + str2double(timeToken{1});
        timeCount = timeCount + 1;
    end
    if ~isempty(rssToken)
        totalRss = totalRss + str2double(rssToken{1});
        rssCount = rssCount + 1;
    end
end

if timeCount == 0
    totalTime = NaN;
end
if rssCount == 0
    totalRss = NaN;
end
end

function curves = sortCurves(curves)
if isempty(curves)
    return;
end

keys = zeros(numel(curves), 2);
for curveI = 1:numel(curves)
    keys(curveI, :) = [levelIndex(curves(curveI).vm), levelIndex(curves(curveI).iion)];
end
[~, order] = sortrows(keys);
curves = curves(order);
end

function label = legendCurveLabel(curve, metricI)
label = sprintf( ...
    '%s | %s/%s | %s | %s', ...
    modelLabel(curve), ...
    formatOrder(curve.ordersAll(metricI)), ...
    formatOrder(curve.ordersFine(metricI)), ...
    formatTime(curve.totalTime), ...
    formatMemory(curve.totalRss) ...
);
end

function label = modelLabel(curve)
label = sprintf('%s-%s', curve.vm, curve.iion);
end

function idx = levelIndex(level)
levels = {'NO', 'p1', 'p2', 'p3'};
idx = find(strcmp(levels, level), 1);
if isempty(idx)
    idx = numel(levels) + 1;
end
end

function color = vmColor(level)
switch level
    case 'NO'
        color = [0.0000 0.4470 0.7410];
    case 'p1'
        color = [0.8500 0.3250 0.0980];
    case 'p2'
        color = [0.4660 0.6740 0.1880];
    case 'p3'
        color = [0.4940 0.1840 0.5560];
    otherwise
        color = [0.2500 0.2500 0.2500];
end
end

function marker = iionMarker(level)
switch level
    case 'NO'
        marker = 'o';
    case 'p1'
        marker = 's';
    case 'p2'
        marker = '^';
    case 'p3'
        marker = 'd';
    otherwise
        marker = 'x';
end
end

function label = alphaLabel(alphaTag)
switch alphaTag
    case '0p0e00'
        label = '0';
    case '1p0em01'
        label = '0.1';
    otherwise
        label = alphaTag;
end
end

function label = dtLabel(dtTag)
switch dtTag
    case '1p0em02'
        label = '1e-2';
    case '1p0em03'
        label = '1e-3';
    otherwise
        label = dtTag;
end
end

function label = meshLabel(meshTag)
switch meshTag
    case 'triangular_Unstr'
        label = 'triangular unstr';
    otherwise
        label = strrep(meshTag, '_', ' ');
end
end

function textValue = formatOrder(value)
if ~isfinite(value)
    textValue = '--';
else
    textValue = sprintf('%.2f', value);
end
end

function textValue = formatTime(seconds)
if ~isfinite(seconds)
    textValue = '--';
elseif seconds < 60
    textValue = sprintf('%.0fs', seconds);
elseif seconds < 3600
    textValue = sprintf('%.0fmin', seconds / 60);
else
    textValue = sprintf('%.1fh', seconds / 3600);
end
end

function textValue = formatMemory(mb)
if ~isfinite(mb)
    textValue = '--';
elseif mb >= 1024
    textValue = sprintf('%.1fGB', mb / 1024);
else
    textValue = sprintf('%.0fMB', mb);
end
end

function occupiedBoxes = addInlineCurveLabel(ax, h, e, label, color, curveIndex, occupiedBoxes)
if nargin < 7
    occupiedBoxes = [];
end
valid = isfinite(h) & isfinite(e) & h > 0 & e > 0;
h = h(valid);
e = e(valid);
if numel(h) < 2
    return;
end

logH = log10(h(:));
logE = log10(e(:));
xMin = min(logH); xMax = max(logH);
yMin = min(logE); yMax = max(logE);
xRange = max(xMax - xMin, eps);
yRange = max(yMax - yMin, eps);

labelFractions = [0.82 0.68 0.54 0.40 0.28 0.92 0.76 0.62 0.48 0.34];
labelFractions = circshift(labelFractions, [0, mod(curveIndex - 1, numel(labelFractions))]);
xOffsets = [0.00 0.05 -0.05 0.09 -0.09 0.13 -0.13 0.18 -0.18];
yOffsets = [0.00 0.045 -0.045 0.085 -0.085 0.125 -0.125 0.165 -0.165];

boxWidth = max(0.055 * xRange, (0.011 * numel(label) + 0.020) * xRange);
boxHeight = 0.050 * yRange;
bestScore = Inf;
bestX = logH(end);
bestY = logE(end);
bestBox = [bestX - boxWidth / 2, bestX + boxWidth / 2, bestY - boxHeight / 2, bestY + boxHeight / 2];

for fracI = 1:numel(labelFractions)
    labelIndex = max(1, min(numel(h), round(labelFractions(fracI) * numel(h))));
    baseX = logH(labelIndex);
    baseY = logE(labelIndex);

    for offsetI = 1:numel(xOffsets)
        candX = baseX + xOffsets(offsetI) * xRange;
        candY = baseY + yOffsets(offsetI) * yRange;
        candX = clampScalar(candX, xMin + 0.035 * xRange, xMax - 0.035 * xRange);
        candY = clampScalar(candY, yMin + 0.050 * yRange, yMax - 0.050 * yRange);
        candBox = [candX - boxWidth / 2, candX + boxWidth / 2, candY - boxHeight / 2, candY + boxHeight / 2];
        score = overlapScore(candBox, occupiedBoxes);
        if score < bestScore
            bestScore = score;
            bestX = candX;
            bestY = candY;
            bestBox = candBox;
        end
        if score == 0
            break;
        end
    end
    if bestScore == 0
        break;
    end
end

text(ax, 10^bestX, 10^bestY, label, ...
    'Color', color, ...
    'FontSize', 5.4, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', color, ...
    'LineWidth', 0.30, ...
    'Margin', 0.7, ...
    'Clipping', 'on' ...
);
occupiedBoxes = [occupiedBoxes; bestBox]; %#ok<AGROW>
end

function value = clampScalar(value, lo, hi)
value = min(max(value, lo), hi);
end

function score = overlapScore(box, boxes)
score = 0.0;
if isempty(boxes)
    return;
end

for boxI = 1:size(boxes, 1)
    xOverlap = max(0.0, min(box(2), boxes(boxI, 2)) - max(box(1), boxes(boxI, 1)));
    yOverlap = max(0.0, min(box(4), boxes(boxI, 4)) - max(box(3), boxes(boxI, 3)));
    score = score + xOverlap * yOverlap;
end
end

function saveFigure(fig, pdfFile, pngFile)
try
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
catch
    print(fig, pdfFile, '-dpdf', '-painters');
end

try
    exportgraphics(fig, pngFile, 'Resolution', 220);
catch
    print(fig, pngFile, '-dpng', '-r220');
end
end
