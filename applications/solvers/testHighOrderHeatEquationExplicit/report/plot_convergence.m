function plot_convergence()
% Plot cell-centred L2 and reconstructed G2 convergence figures for the
% explicit high-order heat-equation MMS example5 runs.
close all; clc;
baseRoot = '/home/pablo/OpenFOAM/pablo-v2312/run/tutorials_electro/testHighOrderHeatEquationExplicit/results/example5';
figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end
files = dir(fullfile(baseRoot, '*', 'hoT_*', 'endT_*_dt_*', 'errors_vs_N_dt_*.dat'));
keys = {};
for i = 1:numel(files)
    fullName = fullfile(files(i).folder, files(i).name);
    parts = split(string(fullName), filesep);
    dimTag = char(parts(find(parts == "example5") + 1));
    dtTok = regexp(files(i).name, 'errors_vs_N_dt_(.+)\.dat', 'tokens', 'once');
    if ~isempty(dtTok); keys{end+1} = [dimTag '|' dtTok{1}]; end %#ok<AGROW>
end
keys = unique(keys);
for k = 1:numel(keys)
    parts = split(keys{k}, '|');
    dimTag = parts{1}; dtTag = parts{2};
    runFiles = dir(fullfile(baseRoot, dimTag, 'hoT_*', ['endT_*_dt_' dtTag], ['errors_vs_N_dt_' dtTag '.dat']));
    plotMetric(runFiles, 'L2', sprintf('heat_explicit_%s_T_L2_dt_%s', dimTag, dtTag), '$||T-T^*||_{L_2}$');
    plotMetric(runFiles, 'G2', sprintf('heat_explicit_%s_T_G2_dt_%s', dimTag, dtTag), '$||T_h-T^*||_{G_2}$');
end
end

function plotMetric(files, metric, outBase, yLabelText)
fig = figure('Color','w','Position',[80 80 900 560]); ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
colors = containers.Map({'NO','p1','p2','p3'}, {[0.1 0.1 0.1], [0.84 0.12 0.12], [0.10 0.32 0.85], [0.72 0.24 0.84]});
markers = containers.Map({'NO','p1','p2','p3'}, {'o','s','^','d'});
allH=[]; allE=[];
for i=1:numel(files)
    T = readtable(fullfile(files(i).folder, files(i).name), 'FileType','text', 'CommentStyle','#', 'Delimiter',' ', 'MultipleDelimsAsOne',true, 'ReadVariableNames',false);
    T.Properties.VariableNames = {'N','L1','L2','Linf','error_T','G1','G2','GInf'};
    level = regexp(files(i).folder, 'hoT_([^/\\]+)', 'tokens', 'once'); if isempty(level); level={'unknown'}; end; level=level{1};
    h = 1./T.N; e = T.(metric); valid = isfinite(h)&isfinite(e)&h>0&e>0; h=h(valid); e=e(valid); [h,ord]=sort(h); e=e(ord);
    c=[0.2 0.2 0.2]; if isKey(colors,level); c=colors(level); end
    m='o'; if isKey(markers,level); m=markers(level); end
    loglog(ax,h,e,'-','Color',c,'Marker',m,'LineWidth',1.2,'MarkerSize',4.5,'DisplayName',['hoT ' level]);
    if ~isempty(h); text(ax,h(end),e(end),['  ' level],'Color',c,'FontSize',8,'BackgroundColor','w','Clipping','on'); end
    allH=[allH; h]; allE=[allE; e]; %#ok<AGROW>
end
if isempty(allH); close(fig); return; end
set(ax,'XDir','reverse'); xlim(ax,[min(allH)/10^0.025 max(allH)*10^0.025]); ylim(ax,10.^[floor(log10(min(allE))) ceil(log10(max(allE)))]);
addSlopes(ax,allH,allE,[2 3 4]); xlabel(ax,'$h=1/N$','Interpreter','latex'); ylabel(ax,yLabelText,'Interpreter','latex'); legend(ax,'Location','eastoutside'); set(ax,'TickLabelInterpreter','latex');
exportgraphics(fig, fullfile(fileparts(mfilename('fullpath')), 'figures', [outBase '.pdf']), 'ContentType','vector');
exportgraphics(fig, fullfile(fileparts(mfilename('fullpath')), 'figures', [outBase '.png']), 'Resolution',300); close(fig);
end

function addSlopes(ax,h,e,orders)
xmin=min(h); xmax=max(h); yy=ylim(ax); ymin=yy(1); ymax=yy(2); hline=logspace(log10(xmin),log10(xmax),100);
for k=1:numel(orders)
    p=orders(k); frac=[0.72 0.52 0.36]; y0=10^(log10(ymin)+frac(k)*(log10(ymax)-log10(ymin))); ref=y0*(hline/xmax).^p;
    loglog(ax,hline,ref,'k:','LineWidth',1.3,'HandleVisibility','off'); tx=xmax/(xmax/xmin)^0.16; ty=y0*(tx/xmax)^p; text(ax,tx,ty,sprintf('$O(h^%d)$',p),'Interpreter','latex','FontSize',8,'BackgroundColor','w','Clipping','on');
end
end
