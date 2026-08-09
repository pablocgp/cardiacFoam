function renderConvergence(dimTag, methodSources, nRangeAll, nRangeFine)
% Shared rendering routine for plot_convergence_2D.m and plot_convergence_3D.m.
%
% Inputs:
%   dimTag        '2D' or '3D' (used in figure-file prefixes and titles).
%   methodSources struct mapping nonlinear-method name -> source dir or
%                 cell-array of source dirs (each containing
%                 hoVm_*_hoStates_* configuration directories that hold
%                 endT_*_dt_<tag>/errors_vs_N_dt_<tag>.dat files). Method
%                 names must match the legacy directory tokens
%                 'Picard', 'diagonalIion', 'JFNK'.
%   nRangeAll     [Nmin Nmax] for the general convergence fit.
%   nRangeFine    [Nmin Nmax] for the asymptotic convergence fit.

figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

methods = fieldnames(methodSources);
methods = orderMethods(methods);
schemes = {'backwardEuler', 'crankNicolson'};
masses = {'consistent', 'lumped'};
fields = {'Vm', 'u1', 'u2'};
fieldText = {'V_m', 'u_1', 'u_2'};
metricNames = {
    {'Vm_L2', 'VmG_L2'};
    {'u1_L2', 'u1G_L2'};
    {'u2_L2', 'u2G_L2'};
};
metricSuffix = {'L2', 'G2'};
metricText = {'cell-centred L_2', 'reconstructed G_2'};

for imethod = 1:numel(methods)
    method = methods{imethod};
    rootList = sourceDirs(methodSources, method);
    if isempty(rootList)
        continue;
    end
    dtTags = findDtTags(rootList, method);

    for ischeme = 1:numel(schemes)
        scheme = schemes{ischeme};
        for imass = 1:numel(masses)
            mass = masses{imass};
            commonAxes = precomputeCommonAxes(rootList, dtTags, method, scheme, mass, metricNames);

            for idt = 1:numel(dtTags)
                dtTag = dtTags{idt};
                files = findErrorFiles(rootList, dtTag, method, scheme, mass);
                if isempty(files)
                    continue;
                end

                for ifield = 1:numel(fields)
                    for imetric = 1:2
                        fig = figure('Color', 'w', 'Position', [80 80 980 640]);
                        ax = axes(fig);
                        hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
                        metric = metricNames{ifield}{imetric};
                        panelTitle = sprintf('%s %s error', fieldText{ifield}, metricText{imetric});
                        plotMetricPanel(ax, files, metric, panelTitle, nRangeAll, nRangeFine);
                        applyCommonAxes(ax, commonAxes{ifield, imetric});
                        addCommonReferenceSlopes(ax, commonAxes{ifield, imetric}, [2 3 4]);
                        title(ax, sprintf('Implicit FDA JfNK: %s, %s, %s, %s mass, %s, dt = %s', ...
                            dimTag, method, schemeShort(scheme), mass, panelTitle, prettyDt(dtTag)), ...
                            'Interpreter', 'tex');

                        outBase = fullfile(figDir, sprintf('implicit_%s_%s_%s_%s_%s_%s_dt_%s', ...
                            dimTag, method, scheme, mass, fields{ifield}, metricSuffix{imetric}, dtTag));
                        exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
                        exportgraphics(fig, [outBase '.png'], 'Resolution', 300);
                        close(fig);
                    end
                end
            end
        end
    end
end
end

function methods = orderMethods(methods)
preferred = {'Picard', 'diagonalIion', 'JFNK'};
keep = preferred(ismember(preferred, methods));
extras = setdiff(methods(:)', preferred, 'stable');
methods = [keep, extras];
end

function dirs = sourceDirs(methodSources, method)
value = methodSources.(method);
if ischar(value) || isstring(value)
    dirs = {char(value)};
elseif iscell(value)
    dirs = cellfun(@char, value, 'UniformOutput', false);
else
    dirs = {};
end
dirs = dirs(cellfun(@(p) ~isempty(p) && exist(p, 'dir'), dirs));
end

function dtTags = findDtTags(rootList, method)
tags = {};
for r = 1:numel(rootList)
    pattern = ['*nonlin_' method '_ode_*'];
    listing = dir(fullfile(rootList{r}, pattern, 'endT_*_dt_*', 'errors_vs_N_dt_*.dat'));
    for i = 1:numel(listing)
        token = regexp(listing(i).folder, 'dt_([^/\\]+)$', 'tokens', 'once');
        if ~isempty(token)
            tags{end+1} = token{1}; %#ok<AGROW>
        end
    end
end
dtTags = unique(tags);
dtTags = sortDtTags(dtTags);
end

function tags = sortDtTags(tags)
if isempty(tags)
    return;
end
vals = zeros(size(tags));
for i = 1:numel(tags)
    vals(i) = str2double(strrep(strrep(tags{i}, 'p', '.'), 'em', 'e-'));
end
[~, order] = sort(vals, 'descend');
tags = tags(order);
end

function files = findErrorFiles(rootList, dtTag, method, scheme, mass)
files = {};
seenConfigs = {};
pattern = ['*scheme_' scheme '_mass_' mass '_nonlin_' method '_ode_*'];
for r = 1:numel(rootList)
    listing = dir(fullfile(rootList{r}, pattern, ['endT_*_dt_' dtTag], ['errors_vs_N_dt_' dtTag '.dat']));
    for i = 1:numel(listing)
        configDir = fileparts(listing(i).folder);
        [~, configName] = fileparts(configDir);
        if any(strcmp(seenConfigs, configName))
            continue;
        end
        seenConfigs{end+1} = configName; %#ok<AGROW>
        files{end+1} = fullfile(listing(i).folder, listing(i).name); %#ok<AGROW>
    end
end
files = sort(files);
end

function commonAxes = precomputeCommonAxes(rootList, dtTags, method, scheme, mass, metricNames)
commonAxes = cell(size(metricNames, 1), 2);
for ifield = 1:size(metricNames, 1)
    for imetric = 1:2
        metric = metricNames{ifield}{imetric};
        allH = [];
        allE = [];
        for idt = 1:numel(dtTags)
            files = findErrorFiles(rootList, dtTags{idt}, method, scheme, mass);
            for i = 1:numel(files)
                T = readErrorFile(files{i});
                if ~ismember(metric, T.Properties.VariableNames)
                    continue;
                end
                h = meshSize(T);
                e = T.(metric);
                valid = isfinite(h) & isfinite(e) & h > 0 & e > 0;
                allH = [allH; h(valid)]; %#ok<AGROW>
                allE = [allE; e(valid)]; %#ok<AGROW>
            end
        end
        commonAxes{ifield, imetric} = niceLogLimits(allH, allE);
    end
end
end

function lims = niceLogLimits(h, e)
valid = isfinite(h) & isfinite(e) & h > 0 & e > 0;
h = h(valid);
e = e(valid);
if isempty(h) || isempty(e)
    lims = [];
    return;
end
xPad = 10^0.025;
xLims = [min(h)/xPad max(h)*xPad];
yLims = 10.^[floor(log10(min(e))) ceil(log10(max(e)))];
if yLims(1) == yLims(2)
    yLims = yLims .* [0.8 1.2];
end
lims = [xLims yLims];
end

function applyCommonAxes(ax, lims)
if isempty(lims)
    return;
end
xlim(ax, lims(1:2));
ylim(ax, lims(3:4));
set(ax, 'XDir', 'reverse');
end

function addCommonReferenceSlopes(ax, lims, orders)
if isempty(lims)
    return;
end
xMin = lims(1); xMax = lims(2); yMin = lims(3); yMax = lims(4);
if any(~isfinite(lims)) || xMin <= 0 || xMax <= 0 || yMin <= 0 || yMax <= 0
    return;
end
hLine = logspace(log10(xMin), log10(xMax), 120);
yStartFractions = [0.72 0.52 0.36];
logYMin = log10(yMin); logYMax = log10(yMax);
for k = 1:numel(orders)
    p = orders(k);
    frac = yStartFractions(min(k, numel(yStartFractions)));
    yStart = 10^(logYMin + frac*(logYMax - logYMin));
    ref = yStart * (hLine / xMax).^p;
    loglog(ax, hLine, ref, 'k:', 'LineWidth', 1.8, 'HandleVisibility', 'off');
    labelX = xMax / (xMax/xMin)^0.14;
    labelY = yStart * (labelX / xMax)^p;
    if labelY > yMin && labelY < yMax
        text(ax, labelX, labelY, sprintf('O(h^%d)', p), ...
            'Interpreter', 'tex', 'FontSize', 7, 'Color', 'k', ...
            'BackgroundColor', 'w', 'Margin', 1, 'Clipping', 'on');
    end
end
end

function plotMetricPanel(ax, files, metric, panelTitle, nRangeAll, nRangeFine)
inlineH = {}; inlineE = {}; inlineLabels = {}; inlineColors = zeros(numel(files), 3);
inlineCount = 0;
legendHandles = gobjects(0); legendLabels = {};
for i = 1:numel(files)
    T = readErrorFile(files{i});
    if ~ismember(metric, T.Properties.VariableNames)
        warning('Metric %s not found in %s', metric, files{i});
        continue;
    end
    h = meshSize(T);
    e = T.(metric);
    valid = isfinite(h) & isfinite(e) & e > 0;
    h = h(valid); e = e(valid);
    [h, order] = sort(h); e = e(order);
    if isempty(h); continue; end
    label = configLabel(files{i});
    [vmLevel, iionLevel] = configLevels(files{i});
    curveColor = vmColor(vmLevel);
    curveMarker = iionMarker(iionLevel);
    [orderAll, orderFine] = convergenceOrders(T, metric, nRangeAll, nRangeFine);
    [totalSeconds, peakRssMB] = accumulatedCost(files{i}, T);
    hp = loglog(ax, h, e, '-', 'Color', curveColor, 'Marker', curveMarker, ...
        'LineWidth', 1.15, 'MarkerSize', 4.5);
    legendHandles(end+1) = hp; %#ok<AGROW>
    legendLabels{end+1} = legendCurveLabel(label, orderAll, orderFine, totalSeconds, peakRssMB); %#ok<AGROW>
    inlineCount = inlineCount + 1;
    inlineH{inlineCount} = h; %#ok<AGROW>
    inlineE{inlineCount} = e; %#ok<AGROW>
    inlineLabels{inlineCount} = shortConfigLabel(files{i}); %#ok<AGROW>
    inlineColors(inlineCount,:) = curveColor;
end
occupiedBoxes = [];
for ilabel = 1:inlineCount
    occupiedBoxes = addInlineCurveLabel(ax, inlineH{ilabel}, inlineE{ilabel}, inlineLabels{ilabel}, inlineColors(ilabel,:), ilabel, occupiedBoxes);
end
xlabel(ax, 'h = 1/N');
ylabel(ax, panelTitle, 'Interpreter', 'tex');
set(ax, 'XScale', 'log', 'YScale', 'log');
if ~isempty(legendHandles)
    lgd = legend(ax, legendHandles, legendLabels, 'Location', 'southwest', 'Interpreter', 'none', 'FontSize', 7.2);
    lgdTitle = title(lgd, 'config | OC all/fine | t_\Sigma | RSS_\Sigma');
    lgdTitle.Interpreter = 'tex';
    if isprop(lgd, 'ItemTokenSize')
        lgd.ItemTokenSize = [10 6];
    end
    lgd.Box = 'on';
    lgd.Color = 'white';
    lgd.EdgeColor = [0.75 0.75 0.75];
end
end

function T = readErrorFile(fileName)
fid = fopen(fileName, 'r');
header = fgetl(fid);
fclose(fid);
header = regexprep(strtrim(header), '^#\s*', '');
names = regexp(header, '\s+', 'split');
varNames = matlab.lang.makeValidName(names);
opts = delimitedTextImportOptions('NumVariables', numel(varNames));
opts.DataLines = [2 Inf];
opts.Delimiter = ' ';
opts.ConsecutiveDelimitersRule = 'join';
opts.LeadingDelimitersRule = 'ignore';
opts.ExtraColumnsRule = 'ignore';
opts.EmptyLineRule = 'read';
opts.VariableNames = varNames;
opts.VariableTypes = repmat({'double'}, 1, numel(varNames));
opts.VariableTypes{1} = 'string';
opts.VariableTypes{2} = 'string';
opts.VariableTypes{5} = 'string';
T = readtable(fileName, opts);
end

function h = meshSize(T)
N = meshN(T);
h = 1 ./ N;
end

function N = meshN(T)
N = double(T.N);
if ismember('dat_file', T.Properties.VariableNames)
    fileNames = string(T.dat_file);
    parsedN = nan(size(N));
    for i = 1:numel(fileNames)
        token = regexp(char(fileNames(i)), '_(\d+)_cells', 'tokens', 'once');
        if ~isempty(token)
            parsedN(i) = str2double(token{1});
        end
    end
    useParsed = isfinite(parsedN) & parsedN > 0;
    N(useParsed) = parsedN(useParsed);
end
end

function [orderAll, orderFine] = convergenceOrders(T, metric, nRangeAll, nRangeFine)
orderAll = NaN;
orderFine = NaN;
if ~ismember(metric, T.Properties.VariableNames)
    return;
end
N = meshN(T);
e = double(T.(metric));
orderAll = fitConvergenceOrder(N, e, nRangeAll(1), nRangeAll(2));
orderFine = fitConvergenceOrder(N, e, nRangeFine(1), nRangeFine(2));
end

function p = fitConvergenceOrder(N, e, nMin, nMax)
p = NaN;
valid = isfinite(N) & isfinite(e) & N >= nMin & N <= nMax & N > 0 & e > 0;
N = N(valid);
e = e(valid);
if numel(N) < 2
    return;
end
h = 1 ./ N;
coeffs = polyfit(log(h), log(e), 1);
p = coeffs(1);
end

function [totalSeconds, peakRssMB] = accumulatedCost(errorFile, T)
totalSeconds = NaN;
peakRssMB = NaN;
if ~ismember('dat_file', T.Properties.VariableNames)
    return;
end
[caseDir, ~, ~] = fileparts(errorFile);
fileNames = string(T.dat_file);
sumTime = 0.0;
sumRss = 0.0;
nTime = 0;
nRss = 0;
for i = 1:numel(fileNames)
    datPath = fullfile(caseDir, char(fileNames(i)));
    if ~exist(datPath, 'file')
        continue;
    end
    txt = fileread(datPath);
    t = scalarFromText(txt, '(?m)^\s*total\s*=\s*([0-9.eE+-]+)\s*$');
    rss = scalarFromText(txt, '(?m)^\s*peakRSS_MB\s*=\s*([0-9.eE+-]+)\s*$');
    if isfinite(t)
        sumTime = sumTime + t;
        nTime = nTime + 1;
    end
    if isfinite(rss)
        sumRss = sumRss + rss;
        nRss = nRss + 1;
    end
end
if nTime > 0
    totalSeconds = sumTime;
end
if nRss > 0
    peakRssMB = sumRss;
end
end

function value = scalarFromText(txt, pattern)
value = NaN;
token = regexp(txt, pattern, 'tokens', 'once');
if isempty(token)
    return;
end
value = str2double(token{1});
end

function label = legendCurveLabel(configLabel, orderAll, orderFine, totalSeconds, peakRssMB)
label = sprintf('%s | %s/%s | %s | %s', ...
    configLabel, formatOrder(orderAll), formatOrder(orderFine), ...
    formatTime(totalSeconds), formatMemory(peakRssMB));
end

function s = formatOrder(value)
if ~isfinite(value)
    s = '--';
else
    s = sprintf('%.2f', value);
end
end

function s = formatTime(seconds)
if ~isfinite(seconds)
    s = '--';
elseif seconds < 60
    s = sprintf('%.0fs', seconds);
elseif seconds < 3600
    s = sprintf('%.0fmin', seconds/60);
else
    s = sprintf('%.0fh', seconds/3600);
end
end

function s = formatMemory(mb)
if ~isfinite(mb)
    s = '--';
elseif mb >= 1024
    s = sprintf('%.0fGB', mb/1024);
else
    s = sprintf('%.0fMB', mb);
end
end

function label = configLabel(fileName)
dtDir = fileparts(fileName); configDir = fileparts(dtDir); [~, configName] = fileparts(configDir);
tokens = regexp(configName, '^hoVm_(.+)_hoStates_(.+)_scheme_(.+)_mass_(.+)_nonlin_(.+)_ode_(.+)$', 'tokens', 'once');
if isempty(tokens)
    label = strrep(configName, 'hoStates', 'hoIion');
else
    label = sprintf('%s-%s', tokens{1}, tokens{2});
end
end

function label = shortConfigLabel(fileName)
label = configLabel(fileName);
end

function [vmLevel, iionLevel] = configLevels(fileName)
dtDir = fileparts(fileName); configDir = fileparts(dtDir); [~, configName] = fileparts(configDir);
tokens = regexp(configName, '^hoVm_(.+)_hoStates_(.+)_scheme_(.+)_mass_(.+)_nonlin_(.+)_ode_(.+)$', 'tokens', 'once');
if isempty(tokens)
    vmLevel = 'unknown'; iionLevel = 'unknown';
else
    vmLevel = tokens{1}; iionLevel = tokens{2};
end
end

function color = vmColor(level)
switch level
    case 'NO'; color = [0.0000 0.4470 0.7410];
    case 'p1'; color = [0.8500 0.3250 0.0980];
    case 'p2'; color = [0.4660 0.6740 0.1880];
    case 'p3'; color = [0.4940 0.1840 0.5560];
    otherwise; color = [0.2500 0.2500 0.2500];
end
end

function marker = iionMarker(level)
switch level
    case 'NO'; marker = 'o';
    case 'p1'; marker = 's';
    case 'p2'; marker = '^';
    case 'p3'; marker = 'd';
    otherwise; marker = 'x';
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
labelFractions = circshift(labelFractions, [0, mod(curveIndex-1, numel(labelFractions))]);
xOffsets = [0.00 0.05 -0.05 0.09 -0.09 0.13 -0.13 0.18 -0.18];
yOffsets = [0.00 0.045 -0.045 0.085 -0.085 0.125 -0.125 0.165 -0.165];

boxWidth = max(0.055*xRange, (0.011*numel(label) + 0.020)*xRange);
boxHeight = 0.050*yRange;
bestScore = Inf;
bestX = logH(end);
bestY = logE(end);
bestBox = [bestX-boxWidth/2, bestX+boxWidth/2, bestY-boxHeight/2, bestY+boxHeight/2];

for ifrac = 1:numel(labelFractions)
    labelIndex = max(1, min(numel(h), round(labelFractions(ifrac)*numel(h))));
    baseX = logH(labelIndex);
    baseY = logE(labelIndex);

    for ioff = 1:numel(xOffsets)
        candX = baseX + xOffsets(ioff)*xRange;
        candY = baseY + yOffsets(ioff)*yRange;
        candX = clampScalar(candX, xMin + 0.035*xRange, xMax - 0.035*xRange);
        candY = clampScalar(candY, yMin + 0.050*yRange, yMax - 0.050*yRange);
        candBox = [candX-boxWidth/2, candX+boxWidth/2, candY-boxHeight/2, candY+boxHeight/2];
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

text(ax, 10^bestX, 10^bestY, label, 'Color', color, 'FontSize', 5.5, 'FontWeight', 'bold', ...
    'Interpreter', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w', 'EdgeColor', color, 'LineWidth', 0.30, 'Margin', 0.7, 'Clipping', 'on');
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
for i = 1:size(boxes, 1)
    xOverlap = max(0.0, min(box(2), boxes(i,2)) - max(box(1), boxes(i,1)));
    yOverlap = max(0.0, min(box(4), boxes(i,4)) - max(box(3), boxes(i,3)));
    score = score + xOverlap*yOverlap;
end
end

function s = prettyDt(tag)
s = strrep(tag, 'p', '.');
s = strrep(s, 'em', 'e-');
end

function s = schemeShort(scheme)
if strcmp(scheme, 'backwardEuler')
    s = 'bE';
elseif strcmp(scheme, 'crankNicolson')
    s = 'cN';
else
    s = scheme;
end
end
