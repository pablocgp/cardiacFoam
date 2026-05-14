
function plot_convergence()
% Plot cell-centred L2 convergence for the implicit manufactured FDA runs.
% This script compares backward Euler and Crank-Nicolson for each mass matrix.

close all; clc;
resultRoot = '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/highOrderManufacturedFDAImplicit/results/example0/2D';
figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

dtTags = {'1p0em03', '1p0em04'};
fields = {'Vm', 'u1', 'u2'};
metrics = containers.Map(fields, {'Vm_L2', 'u1_L2', 'u2_L2'});
pLevels = {'p1', 'p2', 'p3'};
schemes = {'backwardEuler', 'crankNicolson'};
masses = {'consistent', 'lumped'};
lineStyles = containers.Map(schemes, {'-', '-'});
colors = lines(numel(pLevels));
schemeMarkers = containers.Map(schemes, {'o', 's'});

for idt = 1:numel(dtTags)
    dtTag = dtTags{idt};
    for ifield = 1:numel(fields)
        fld = fields{ifield};
        metric = metrics(fld);
        fig = figure('Color', 'w', 'Position', [60 80 1450 650]);
        tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        for imass = 1:numel(masses)
            mass = masses{imass};
            ax = nexttile; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
            allH = []; allE = [];

            for ip = 1:numel(pLevels)
                p = pLevels{ip};
                for is = 1:numel(schemes)
                    scheme = schemes{is};
                    cfg = sprintf('hoVm_%s_hoStates_%s_scheme_%s_mass_%s', p, p, scheme, mass);
                    fileName = fullfile(resultRoot, cfg, ['endT_0.05_dt_' dtTag], ['errors_vs_N_dt_' dtTag '.dat']);
                    if ~exist(fileName, 'file')
                        warning('Missing %s', fileName);
                        continue;
                    end
                    T = readErrorFile(fileName);
                    if ~ismember(metric, T.Properties.VariableNames)
                        warning('Metric %s not found in %s', metric, fileName);
                        continue;
                    end
                    h = 1 ./ T.N;
                    e = T.(metric);
                    valid = isfinite(h) & isfinite(e) & e > 0;
                    h = h(valid); e = e(valid);
                    [h, order] = sort(h);
                    e = e(order);
                    loglog(ax, h, e, lineStyles(scheme), ...
                        'Color', colors(ip,:), ...
                        'Marker', schemeMarkers(scheme), ...
                        'LineWidth', 1.35, ...
                        'MarkerSize', 4.5, ...
                        'DisplayName', sprintf('%s %s', p, scheme));
                    allH = [allH; h(:)]; %#ok<AGROW>
                    allE = [allE; e(:)]; %#ok<AGROW>
                end
            end

            addReferenceSlopes(ax, allH, allE, [2 3 4]);
            xlabel(ax, 'h = 1/N');
            ylabel(ax, sprintf('%s cell-centred L2 error', fld));
            title(ax, sprintf('%s mass', mass), 'Interpreter', 'none');
            set(ax, 'XScale', 'log', 'YScale', 'log');
            legend(ax, 'Location', 'best', 'Interpreter', 'none', 'FontSize', 7);
        end

        sgtitle(fig, sprintf('Implicit manufactured FDA: %s, dt = %s', fld, prettyDt(dtTag)), 'Interpreter', 'none');
        outBase = fullfile(figDir, sprintf('implicit_%s_dt_%s', fld, dtTag));
        exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
        exportgraphics(fig, [outBase '.png'], 'Resolution', 300);
    end
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
