function experimentHuman(cfg)
% EXPERIMENTHUMAN  Part 18 / Fig. 7-8: synthetic walking-human point
% cloud placed at the paper's three lateral positions (y = 1.1, 0,
% -0.8 m, z = 1 m standoff), compared across:
%   (a) Planar-scanning SISO-SAR + classical BPA
%   (b) Proposed curved SISO-SAR + classical BPA
%   (c) Proposed curved SISO-SAR + proposed measurement-fusion algorithm
%
% AXIS CONVENTION (see targets/createHumanPointCloud.m): x = height,
% y = lateral position, z = depth. For each of the three positions we
% render an (x,z) slice through that position's y -- i.e. a vertical
% profile of the person at their standoff depth. The three positions
% are labeled "front"/"side"/"back" only loosely, following the
% paper's own framing that the array's wide FoV lets different lateral
% placements present different aspects of the body to the array; we do
% NOT claim to reproduce the paper's exact (unstated) body-relative
% imaging-plane convention. This is our own reasonable, self-consistent
% simplification, documented here and in the README.

if nargin < 1 || isempty(cfg)
    cfg = config();
end

fprintf('\n=== Experiment 5: Human Point Cloud ===\n');

curvedAp = generateCurvedArray(cfg);
planarAp = generatePlanarArray(cfg);
[~, k] = generateFrequencyGrid(cfg);

yOffsets = [1.1, 0, -0.8];
zDepth = 1.0;
viewNames = {'y=1.1m ("front")', 'y=0m ("side")', 'y=-0.8m ("back")'};

res = cfg.imaging.voxelRes;
xAxis = -1.0:res:0.95;                 % height axis, covers the whole body
zAxis = zDepth + (-0.3:res:0.3);       % depth window around the body's thickness

methodNames = {'Planar SISO-SAR + BPA', 'Curved SISO-SAR + BPA', 'Curved SISO-SAR + Proposed'};
imgCell = cell(3,3);   % {method, position}

for p = 1:3
    human = createHumanPointCloud(yOffsets(p), zDepth, 1);
    fprintf('Position %d/%d (%s): %d scatterers\n', p, 3, viewNames{p}, size(human.pos,1));

    [X, Z] = meshgrid(xAxis, zAxis);
    voxelPos = [X(:), yOffsets(p)*ones(numel(X),1), Z(:)];

    echoPlanar = generateHumanPointCloudEcho(planarAp.pos, human, k, planarAp.beamDir, cfg.HPBW_deg);
    echoCurved = generateHumanPointCloudEcho(curvedAp.pos, human, k, curvedAp.beamDir, cfg.HPBW_deg);
    echoPlanar = maybeAddNoise(echoPlanar, cfg);   % OUR EXTENSION: no-op unless cfg.noise.enable
    echoCurved = maybeAddNoise(echoCurved, cfg);

    imgA = BPA(echoPlanar, planarAp.pos, k, voxelPos, cfg);
    imgB = BPA(echoCurved, curvedAp.pos, k, voxelPos, cfg);
    imgC = measurementFusionReconstruction(curvedAp, echoCurved, k, voxelPos, cfg);

    imgCell{1,p} = reshape(toDB(imgA, cfg.imaging.dB_range), size(X));
    imgCell{2,p} = reshape(toDB(imgB, cfg.imaging.dB_range), size(X));
    imgCell{3,p} = reshape(toDB(imgC, cfg.imaging.dB_range), size(X));
end

% ---- One figure per method, 3 columns (positions) ----
fnames = {'human_planar.png', 'human_curved_BPA.png', 'human_curved_proposed.png'};
for m = 1:3
    plotComparison(imgCell(m,:), viewNames, xAxis, zAxis, 'x / m (height)', 'z / m (depth)', ...
        cfg.imaging.dB_range, methodNames{m}, fullfile(cfg.resultsDir, fnames{m}));
end

% ---- One combined 3x3 figure (rows = method, cols = position) ----
fig = figure('Visible', 'on', 'Position', [50 50 1200 1150]);
drawnow;
tl = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for m = 1:3
    for p = 1:3
        nexttile(tl, (m-1)*3+p);
        imagesc(xAxis, zAxis, imgCell{m,p});
        axis xy equal tight; colormap('turbo'); caxis([-cfg.imaging.dB_range, 0]); colorbar;
        xlabel('x/m'); ylabel('z/m');
        title(sprintf('%s - %s', methodNames{m}, viewNames{p}), 'Interpreter', 'none', 'FontSize', 8);
    end
end
title(tl, 'Walking-Human Point Cloud: System/Algorithm Comparison');
exportgraphics(fig, fullfile(cfg.resultsDir, 'human_full_comparison.png'), 'Resolution', 150);

save(fullfile(cfg.resultsDir, 'human_data.mat'), 'imgCell', 'xAxis', 'zAxis', 'yOffsets', 'zDepth', 'methodNames', 'viewNames');

fprintf('Saved: human_planar.png, human_curved_BPA.png, human_curved_proposed.png, human_full_comparison.png, human_data.mat\n');
end
