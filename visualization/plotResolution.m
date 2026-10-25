function fig = plotResolution(angleDeg, seriesCell, legendCell, titleStr, ylabelStr, savePath)
% PLOTRESOLUTION  Resolution-vs-elevation-angle comparison across
% systems/algorithms (Fig. 5 style). Each series is in meters; plotted
% in cm to match the paper's axis units.

fig = figure('Visible', 'on');
drawnow;
hold on; grid on; box on;
markers = {'-o', '-s', '-^', '-d'};
for i = 1:numel(seriesCell)
    m = markers{mod(i-1, numel(markers)) + 1};
    plot(angleDeg, seriesCell{i} * 100, m, 'LineWidth', 1.4, 'MarkerSize', 5);
end
xlabel('angle / deg');
ylabel(ylabelStr);
title(titleStr, 'Interpreter', 'none');
legend(legendCell, 'Location', 'best');
set(gca, 'FontSize', 10);

if nargin >= 6 && ~isempty(savePath)
    exportgraphics(fig, savePath, 'Resolution', 150);
end
end
