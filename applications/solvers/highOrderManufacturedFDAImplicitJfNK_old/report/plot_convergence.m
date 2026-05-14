function plot_convergence()
% Plot individual cell-centred L2 and reconstructed G2 convergence figures
% for the implicit manufactured FDA runs.
%
% One figure file is generated for each dimension, time step, time scheme,
% mass matrix, variable, and error metric. The legend is placed outside the
% axes so it does not hide the convergence curves.
%
% NOTE: result directories still use the legacy token "hoStates". In the
% current solver this token denotes the high-order quadrature level used for
% Iion and for the local ODE state values stored at the Iion integration
% points.

close all; clc;
baseRoot = '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicit/results/example0';
figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

dimTags = listSubdirs(baseRoot);
dtTags = findDtTags(baseRoot);
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

for idim = 1:numel(dimTags)
    dimTag = dimTags{idim};
    resultRoot = fullfile(baseRoot, dimTag);

    for ischeme = 1:numel(schemes)
        scheme = schemes{ischeme};
        for imass = 1:numel(masses)
            mass = masses{imass};
            commonAxes = precomputeCommonAxes(resultRoot, dtTags, scheme, mass, metricNames);

            for idt = 1:numel(dtTags)
                dtTag = dtTags{idt};
                files = findErrorFiles(resultRoot, dtTag, scheme, mass);
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
                        plotMetricPanel(ax, files, metric, panelTitle);
                        applyCommonAxes(ax, commonAxes{ifield, imetric});
                        addCommonReferenceSlopes(ax, commonAxes{ifield, imetric}, [2 3 4]);
                        title(ax, sprintf('Implicit FDA: %s, %s, %s mass, %s, dt = %s', ...
                            dimTag, schemeShort(scheme), mass, panelTitle, prettyDt(dtTag)), ...
                            'Interpreter', 'tex');

                        outBase = fullfile(figDir, sprintf('implicit_%s_%s_%s_%s_%s_dt_%s', ...
                            dimTag, scheme, mass, fields{ifield}, metricSuffix{imetric}, dtTag));
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

function dirs = listSubdirs(rootDir)
listing = dir(rootDir);
dirs = {listing([listing.isdir]).name};
dirs = dirs(~ismember(dirs, {'.', '..'}));
dirs = sort(dirs);
end

function dtTags = findDtTags(baseRoot)
listing = dir(fullfile(baseRoot, '*', '*', 'endT_*_dt_*', 'errors_vs_N_dt_*.dat'));
tags = {};
for i = 1:numel(listing)
    token = regexp(listing(i).folder, 'dt_([^/\\]+)$', 'tokens', 'once');
    if ~isempty(token)
        tags{end+1} = token{1}; %#ok<AGROW>
    end
end
dtTags = unique(tags);
dtTags = sortDtTags(dtTags);
end

function tags = sortDtTags(tags)
vals = zeros(size(tags));
for i = 1:numel(tags)
    vals(i) = str2double(strrep(strrep(tags{i}, 'p', '.'), 'em', 'e-'));
end
[~, order] = sort(vals, 'descend');
tags = tags(order);
end

function files = findErrorFiles(resultRoot, dtTag, scheme, mass)
pattern = ['*scheme_' scheme '_mass_' mass];
listing = dir(fullfile(resultRoot, pattern, ['endT_*_dt_' dtTag], ['errors_vs_N_dt_' dtTag '.dat']));
files = arrayfun(@(d) fullfile(d.folder, d.name), listing, 'UniformOutput', false);
files = sort(files);
end

function commonAxes = precomputeCommonAxes(resultRoot, dtTags, scheme, mass, metricNames)
commonAxes = cell(size(metricNames, 1), 2);
for ifield = 1:size(metricNames, 1)
    for imetric = 1:2
        metric = metricNames{ifield}{imetric};
        allH = [];
        allE = [];
        for idt = 1:numel(dtTags)
            files = findErrorFiles(resultRoot, dtTags{idt}, scheme, mass);
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

function plotMetricPanel(ax, files, metric, panelTitle)
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
    hp = loglog(ax, h, e, '-', 'Color', curveColor, 'Marker', curveMarker, ...
        'LineWidth', 1.15, 'MarkerSize', 4.5, 'DisplayName', label);
    legendHandles(end+1) = hp; %#ok<AGROW>
    legendLabels{end+1} = label; %#ok<AGROW>
    inlineCount = inlineCount + 1;
    inlineH{inlineCount} = h; %#ok<AGROW>
    inlineE{inlineCount} = e; %#ok<AGROW>
    inlineLabels{inlineCount} = shortConfigLabel(files{i}); %#ok<AGROW>
    inlineColors(inlineCount,:) = curveColor;
end
for ilabel = 1:inlineCount
    addInlineCurveLabel(ax, inlineH{ilabel}, inlineE{ilabel}, inlineLabels{ilabel}, inlineColors(ilabel,:), ilabel);
end
xlabel(ax, 'h = 1/N');
ylabel(ax, panelTitle, 'Interpreter', 'tex');
set(ax, 'XScale', 'log', 'YScale', 'log');
if ~isempty(legendHandles)
    lgd = legend(ax, legendHandles, legendLabels, 'Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 7);
    lgd.Box = 'off';
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
N = T.N;
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
h = 1 ./ N;
end

function label = configLabel(fileName)
dtDir = fileparts(fileName); configDir = fileparts(dtDir); [~, configName] = fileparts(configDir);
tokens = regexp(configName, '^hoVm_(.+)_hoStates_(.+)_scheme_(.+)_mass_(.+)$', 'tokens', 'once');
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
tokens = regexp(configName, '^hoVm_(.+)_hoStates_(.+)_scheme_(.+)_mass_(.+)$', 'tokens', 'once');
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

function addInlineCurveLabel(ax, h, e, label, color, curveIndex)
if numel(h) < 2 || isempty(e); return; end
labelFractions = [0.78 0.64 0.50 0.36 0.86 0.72 0.58 0.44 0.30 0.92 0.68 0.54];
labelFraction = labelFractions(mod(curveIndex-1, numel(labelFractions)) + 1);
labelIndex = max(1, min(numel(h), round(labelFraction*numel(h))));
x = h(labelIndex); y = e(labelIndex);
if ~isfinite(x) || ~isfinite(y) || x <= 0 || y <= 0; return; end
xOffsetCycle = [1.18 0.86 1.34 0.74 1.55 0.64 1.78 0.56 1.08 0.92 1.24 0.80];
yOffsetCycle = [1.18 0.82 1.38 0.68 1.62 0.56 1.92 0.48 1.08 0.92 1.28 0.76];
x = x*xOffsetCycle(mod(curveIndex-1, numel(xOffsetCycle)) + 1);
y = y*yOffsetCycle(mod(curveIndex-1, numel(yOffsetCycle)) + 1);
xMin = min(h); xMax = max(h);
x = min(max(x, 1.03*xMin), 0.97*xMax);
text(ax, x, y, label, 'Color', color, 'FontSize', 6, 'FontWeight', 'bold', ...
    'Interpreter', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w', 'EdgeColor', color, 'LineWidth', 0.35, 'Margin', 1, 'Clipping', 'on');
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
