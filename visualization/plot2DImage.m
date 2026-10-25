function fig = plot2DImage(xAxis, yAxis, imgDB, titleStr, xlabelStr, ylabelStr, dynRangeDB, savePath)
% PLOT2DIMAGE  Single normalized-dB 2-D reconstruction image.
%   imgDB : length(yAxis) x length(xAxis) matrix (rows = y, cols = x --
%           i.e. the layout produced by meshgrid(xAxis,yAxis)).

fig = figure('Visible', 'on');
drawnow;
imagesc(xAxis, yAxis, imgDB);
axis xy equal tight;
colormap('turbo');
caxis([-dynRangeDB, 0]);
cb = colorbar; cb.Label.String = 'dB';
xlabel(xlabelStr); ylabel(ylabelStr);
title(titleStr, 'Interpreter', 'none');
set(gca, 'FontSize', 10);

if nargin >= 8 && ~isempty(savePath)
    exportgraphics(fig, savePath, 'Resolution', 150);
end
end
