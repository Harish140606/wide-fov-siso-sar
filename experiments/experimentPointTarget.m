function experimentPointTarget(cfg)
% EXPERIMENTPOINTTARGET  Part 15: single point-target 3-D reconstruction
% with classical BPA and the proposed measurement-fusion algorithm,
% both on the curved-array system. Verifies the reconstructed peak is
% near the true target location (Part 22 validation).

if nargin < 1 || isempty(cfg)
    cfg = config();
end

fprintf('\n=== Experiment 2: Single Point Target ===\n');

ap = generateCurvedArray(cfg);
[~, k] = generateFrequencyGrid(cfg);

tgt = createPointTarget(cfg.target.r, cfg.target.theta_deg, cfg.target.phi_deg);
fprintf('Target: r=%.2fm theta=%.1fdeg phi=%.1fdeg -> pos=[%.3f %.3f %.3f]\n', ...
    cfg.target.r, cfg.target.theta_deg, cfg.target.phi_deg, tgt.pos);

echo = generatePointTargetEcho(ap.pos, tgt.pos, tgt.amp, k, ap.beamDir, cfg.HPBW_deg);   % Eq. (1)
echo = maybeAddNoise(echo, cfg);   % OUR EXTENSION: no-op unless cfg.noise.enable

hw  = cfg.imaging.windowHalfWidth;
res = cfg.imaging.voxelRes;
xAxis = tgt.pos(1) + (-hw:res:hw);
yAxis = tgt.pos(2) + (-hw:res:hw);
zAxis = tgt.pos(3) + (-hw:res:hw);
[X, Y, Z] = meshgrid(xAxis, yAxis, zAxis);
voxelPos = [X(:), Y(:), Z(:)];
fprintf('3-D imaging grid: %d x %d x %d = %d voxels\n', numel(yAxis), numel(xAxis), numel(zAxis), numel(X));

tic;
imgBPA = BPA(echo, ap.pos, k, voxelPos, cfg);                  % Eq. (8)
tBPA = toc;

tic;
[imgMF, diagMF] = measurementFusionReconstruction(ap, echo, k, voxelPos, cfg);  % Eq. (13)
tMF = toc;

fprintf('Runtime: BPA=%.2fs, Proposed=%.2fs\n', tBPA, tMF);

volBPA = reshape(imgBPA, size(X));
volMF  = reshape(imgMF, size(X));

[~, pkI]  = max(abs(volBPA(:)));
[~, pkI2] = max(abs(volMF(:)));
peakBPA = [X(pkI) Y(pkI) Z(pkI)];
peakMF  = [X(pkI2) Y(pkI2) Z(pkI2)];
errBPA = norm(peakBPA - tgt.pos);
errMF  = norm(peakMF - tgt.pos);
fprintf('Peak location error: BPA=%.4fm, Proposed=%.4fm (voxel size=%.4fm)\n', errBPA, errMF, res);

volBPAdB = reshape(toDB(volBPA(:), cfg.imaging.dB_range), size(X));
volMFdB  = reshape(toDB(volMF(:),  cfg.imaging.dB_range), size(X));

plot3DImage(xAxis, yAxis, zAxis, volBPAdB, -6, ...
    'Point Target - Classical BPA (-6 dB isosurface)', ...
    fullfile(cfg.resultsDir, 'point_target_BPA.png'));
plot3DImage(xAxis, yAxis, zAxis, volMFdB, -6, ...
    'Point Target - Proposed Measurement-Fusion (-6 dB isosurface)', ...
    fullfile(cfg.resultsDir, 'point_target_proposed.png'));

[~, zi] = min(abs(zAxis - tgt.pos(3)));
plotComparison({volBPAdB(:,:,zi), volMFdB(:,:,zi)}, {'Classical BPA', 'Proposed'}, ...
    xAxis, yAxis, 'x / m', 'y / m', cfg.imaging.dB_range, ...
    sprintf('Point-Target Slice at z = %.3f m', zAxis(zi)), ...
    fullfile(cfg.resultsDir, 'point_target_slice_comparison.png'));

save(fullfile(cfg.resultsDir, 'point_target_data.mat'), ...
    'tgt', 'xAxis', 'yAxis', 'zAxis', 'volBPA', 'volMF', 'errBPA', 'errMF', 'diagMF');

fprintf('Saved: point_target_BPA.png, point_target_proposed.png, point_target_slice_comparison.png, point_target_data.mat\n');
end
