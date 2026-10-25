function fig = plotComparison(imgCellDB, titleCell, xAxis, yAxis, xlabelStr, ylabelStr, dynRangeDB, overallTitle, savePath)
% PLOTCOMPARISON  Side-by-side normalized-dB images for N reconstruction
% methods/systems (used for the BPA-vs-proposed and
% planar-vs-curved comparisons, Fig. 6 / Fig. 8 style).

n = numel(imgCellDB);
fig = figure('Visible', 'on', 'Position', [50 50 420*n 420]);
drawnow;
tl = tiledlayout(fig, 1, n, 'Padding', 'compact', 'TileSpacing', 'compact');
for idx = 1:n
    nexttile(tl);
    imagesc(xAxis, yAxis, imgCellDB{idx});
    axis xy equal tight;
    colormap('turbo');
    caxis([-dynRangeDB, 0]);
    colorbar;
    xlabel(xlabelStr); ylabel(ylabelStr);
    title(titleCell{idx}, 'Interpreter', 'none');
end
title(tl, overallTitle, 'Interpreter', 'none');

if nargin >= 9 && ~isempty(savePath)
    exportgraphics(fig, savePath, 'Resolution', 150);
end
end
