function fig = plotSamplingPoints(pos, titleStr, savePath)
% PLOTSAMPLINGPOINTS  3-D scatter of aperture sampling points
% (reproduces the visual intent of Fig. 3(a) / Fig. 4).

fig = figure('Visible', 'on');
drawnow;
scatter3(pos(:,1), pos(:,2), pos(:,3), 8, pos(:,3), 'filled');
axis equal; grid on;
xlabel('x / m'); ylabel('y / m'); zlabel('z / m');
title(titleStr, 'Interpreter', 'none');
view(40, 25);
colormap('turbo');
cb = colorbar; cb.Label.String = 'z / m';

if nargin >= 3 && ~isempty(savePath)
    exportgraphics(fig, savePath, 'Resolution', 150);
end
end
