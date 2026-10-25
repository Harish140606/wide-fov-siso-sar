function fig = plot3DImage(xAxis, yAxis, zAxis, volDB, isoLevelDB, titleStr, savePath)
% PLOT3DIMAGE  Isosurface rendering of a 3-D reconstructed field.
%   volDB : length(yAxis) x length(xAxis) x length(zAxis), matching
%           meshgrid(xAxis,yAxis,zAxis) layout.
%   isoLevelDB : dB level of the isosurface to draw (e.g. -6)

fig = figure('Visible', 'on');
drawnow;
[X, Y, Z] = meshgrid(xAxis, yAxis, zAxis);

if all(volDB(:) <= isoLevelDB)
    % Degenerate case: nothing reaches the requested level. Fall back
    % to plotting the single brightest voxel so the figure is not empty.
    [~, li] = max(volDB(:));
    [iy, ix, iz] = ind2sub(size(volDB), li);
    scatter3(X(iy,ix,iz), Y(iy,ix,iz), Z(iy,ix,iz), 80, 'filled', 'r');
else
    p = patch(isosurface(X, Y, Z, volDB, isoLevelDB));
    isonormals(X, Y, Z, volDB, p);
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.85;
    camlight; lighting gouraud;
end

axis equal; grid on; view(35, 20);
xlabel('x / m'); ylabel('y / m'); zlabel('z / m');
title(titleStr, 'Interpreter', 'none');

if nargin >= 7 && ~isempty(savePath)
    exportgraphics(fig, savePath, 'Resolution', 150);
end
end
