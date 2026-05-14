
function plot_convergence()
% Plot cell-centred L2 convergence for the explicit manufactured FDA runs.
% This script only reads existing .dat files and writes figures/ outputs.

close all; clc;
resultRoot = '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAExplicit/results/example0/2D';
figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

dtTags = {'1p0em03', '1p0em04', '1p0em05'};
fields = {'Vm', 'u1', 'u2'};
metrics = containers.Map(fields, {'Vm_L2', 'u1_L2', 'u2_L2'});

for idt = 1:numel(dtTags)
    dtTag = dtTags{idt};
    files = findErrorFiles(resultRoot, dtTag);
    if isempty(files)
        warning('No errors_vs_N file found for dt tag %s', dtTag);
        continue;
    end

    for ifield = 1:numel(fields)
        fld = fields{ifield};
        metric = metrics(fld);
        fig = figure('Color', 'w', 'Position', [80 80 1320 780]);
        ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        colors = lines(max(numel(files), 1));
        markers = {'o','s','^','d','v','>','<','p','h','x','+'};
        allH = []; allE = [];

        for i = 1:numel(files)
            T = readErrorFile(files{i});
            if ~ismember(metric, T.Properties.VariableNames)
                warning('Metric %s not found in %s', metric, files{i});
                continue;
            end
            h = 1 ./ T.N;
            e = T.(metric);
            valid = isfinite(h) & isfinite(e) & e > 0;
            h = h(valid); e = e(valid);
            [h, order] = sort(h);
            e = e(order);
            label = configLabel(files{i});
            loglog(ax, h, e, '-', ...
                'Color', colors(i,:), ...
                'Marker', markers{mod(i-1, numel(markers))+1}, ...
                'LineWidth', 1.15, ...
                'MarkerSize', 4.5, ...
                'DisplayName', label);
            allH = [allH; h(:)]; %#ok<AGROW>
            allE = [allE; e(:)]; %#ok<AGROW>
        end

        addReferenceSlopes(ax, allH, allE, [2 3 4]);
        xlabel(ax, 'h = 1/N');
        ylabel(ax, sprintf('%s cell-centred L2 error', fld));
        title(ax, sprintf('Explicit manufactured FDA: %s, dt = %s', fld, prettyDt(dtTag)), 'Interpreter', 'none');
        legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 7);
        set(ax, 'XScale', 'log', 'YScale', 'log');

        outBase = fullfile(figDir, sprintf('explicit_%s_dt_%s', fld, dtTag));
        exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
        exportgraphics(fig, [outBase '.png'], 'Resolution', 300);
    end
end
end

function files = findErrorFiles(resultRoot, dtTag)
listing = dir(fullfile(resultRoot, '*', ['endT_*_dt_' dtTag], ['errors_vs_N_dt_' dtTag '.dat']));
files = arrayfun(@(d) fullfile(d.folder, d.name), listing, 'UniformOutput', false);
files = sort(files);
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
opts.VariableTypes{1} = 'string'; % dim
opts.VariableTypes{2} = 'string'; % config
opts.VariableTypes{5} = 'string'; % dat_file
T = readtable(fileName, opts);
end

function label = configLabel(fileName)
dtDir = fileparts(fileName);
configDir = fileparts(dtDir);
[~, label] = fileparts(configDir);
end

function addReferenceSlopes(ax, h, e, orders)
valid = isfinite(h) & isfinite(e) & h > 0 & e > 0;
h = h(valid); e = e(valid);
if isempty(h) || isempty(e); return; end
hLine = logspace(log10(min(h)), log10(max(h)), 120);
h0 = exp(mean(log(h)));
e0 = exp(mean(log(e)));
offsets = logspace(-0.35, 0.35, numel(orders));
for k = 1:numel(orders)
    p = orders(k);
    ref = e0 * offsets(k) * (hLine / h0).^p;
    loglog(ax, hLine, ref, 'k:', 'LineWidth', 1.8, 'DisplayName', sprintf('O(h^%d)', p));
end
end

function s = prettyDt(tag)
s = strrep(tag, 'p', '.');
s = strrep(s, 'em', 'e-');
end
