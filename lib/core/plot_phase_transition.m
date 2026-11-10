function fig = plot_phase_transition(results, threshold, out_base, figcfg)
%PLOT_PHASE_TRANSITION Recovery rate over (alpha, m): one row per dimension d, one column per method.

if nargin < 2 || isempty(threshold), threshold = 1; end
if nargin < 3, out_base = ''; end
if nargin < 4, cfg = setup(); figcfg = cfg.fig; end
if ~istable(results)
    if isfolder(results), results = fullfile(results, 'per_trial.csv'); end
    results = readtable(results, 'TextType', 'char', 'Delimiter', ',');
end
T = results;
T.alpha = round(T.alpha, 6);
order = {'SREDG', 'PREEEDM', 'RGNMR'};
methods = [order(ismember(order, T.method)), setdiff(unique(T.method, 'stable'), order, 'stable')'];
dims = unique(T.d)'; ms = unique(T.m); as = unique(T.alpha);
nr = numel(dims); nc = numel(methods);

fig = figure('Color', 'w', 'Position', [100 100 figcfg.panel_size(1) * nc + 90, figcfg.panel_size(2) * nr + 40]);
tl = tiledlayout(nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
k = 0;
for i = 1:nr
    for j = 1:nc
        k = k + 1;
        ax = nexttile(tl);
        R = T(T.d == dims(i) & strcmp(T.method, methods{j}), :);
        Z = nan(numel(ms), numel(as));                  % rows m, columns alpha
        for a = 1:numel(ms)
            for b = 1:numel(as)
                sel = R.m == ms(a) & R.alpha == as(b);
                if any(sel), Z(a, b) = mean(R.rmse(sel) < threshold); end
            end
        end
        imagesc(ax, 100 * as, ms, Z, 'AlphaData', ~isnan(Z), [0 1]);
        set(ax, 'YDir', 'normal', 'Color', figcfg.nodata_color, 'TickDir', 'in', 'Box', 'on', ...
            'FontSize', figcfg.font_size);
        colormap(ax, feval(figcfg.colormap, 256));
        title(ax, sprintf('(%c) %s, r = %d', 'a' + k - 1, methods{j}, dims(i)), 'FontWeight', 'normal');
        if j == 1, ylabel(ax, 'Number of anchors, m'); else, ax.YTickLabel = []; end
        if i == nr, xlabel(ax, 'Outliers, \alpha (%)'); else, ax.XTickLabel = []; end
    end
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = sprintf('Recovery rate (RMSE < %g)', threshold);
cb.FontSize = figcfg.font_size;

if ~isempty(out_base)
    for f = figcfg.formats
        if strcmp(f{1}, 'png')
            exportgraphics(fig, [out_base '.png'], 'Resolution', figcfg.resolution);
        else
            exportgraphics(fig, [out_base '.' f{1}], 'ContentType', 'vector');
        end
    end
end
end
