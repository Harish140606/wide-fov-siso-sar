function experimentPointArray(cfg)
% EXPERIMENTPOINTARRAY  Part 17 / Fig. 6: planar point-array target at
% z = 1 m, imaged by the curved SISO-SAR system with (a) classical BPA
% and (b) the proposed measurement-fusion algorithm. Demonstrates the
% central circular artifacts / peripheral detail loss / sidelobes
% produced by BPA on the non-uniform aperture, and their mitigation by
% the proposed algorithm.

if nargin < 1 || isempty(cfg)
    cfg = config();
end

fprintf('\n=== Experiment 4: Planar Point Array (z = 1 m) ===\n');

ap = generateCurvedArray(cfg);
[~, k] = generateFrequencyGrid(cfg);

z0 = 1.0;
halfWidth = 0.4;
spacing   = 0.1;
targets = createPointArray(z0, halfWidth, spacing);
fprintf('Point-array target: %d scatterers on a %.1fx%.1fm grid at z=%.1fm\n', ...
    size(targets.pos,1), 2*halfWidth, 2*halfWidth, z0);

echo = generateSceneEcho(ap.pos, targets, k, ap.beamDir, cfg.HPBW_deg);
echo = maybeAddNoise(echo, cfg);   % OUR EXTENSION: no-op unless cfg.noise.enable

hw  = 0.5;
res = cfg.imaging.voxelRes;
xAxis = -hw:res:hw;
yAxis = -hw:res:hw;
[X, Y] = meshgrid(xAxis, yAxis);
voxelPos = [X(:), Y(:), z0*ones(numel(X), 1)];
fprintf('Imaging grid: %d x %d = %d pixels\n', numel(yAxis), numel(xAxis), numel(X));

tic; imgBPA = BPA(echo, ap.pos, k, voxelPos, cfg); tBPA = toc;
tic; [imgMF, diagMF] = measurementFusionReconstruction(ap, echo, k, voxelPos, cfg); tMF = toc;
fprintf('Runtime: BPA=%.2fs, Proposed=%.2fs\n', tBPA, tMF);

volBPAdB = reshape(toDB(imgBPA, cfg.imaging.dB_range), size(X));
volMFdB  = reshape(toDB(imgMF,  cfg.imaging.dB_range), size(X));

plot2DImage(xAxis, yAxis, volBPAdB, 'Curved SISO-SAR + Classical BPA (z = 1 m)', ...
    'x / m', 'y / m', cfg.imaging.dB_range, fullfile(cfg.resultsDir, 'point_array_BPA.png'));
plot2DImage(xAxis, yAxis, volMFdB, 'Curved SISO-SAR + Proposed Algorithm (z = 1 m)', ...
    'x / m', 'y / m', cfg.imaging.dB_range, fullfile(cfg.resultsDir, 'point_array_proposed.png'));
plotComparison({volBPAdB, volMFdB}, {'Classical BPA', 'Proposed Measurement-Fusion'}, ...
    xAxis, yAxis, 'x / m', 'y / m', cfg.imaging.dB_range, ...
    'Planar Point-Array Imaging at z = 1 m', fullfile(cfg.resultsDir, 'point_array_comparison.png'));

% ---- Metrics ----
magBPA = reshape(abs(imgBPA), size(X));
magMF  = reshape(abs(imgMF),  size(X));

trueMask = false(size(X));
tol = res * 1.5;
for t = 1:size(targets.pos, 1)
    d = sqrt((X - targets.pos(t,1)).^2 + (Y - targets.pos(t,2)).^2);
    trueMask = trueMask | (d < tol);
end

artifactBPA = computeArtifactEnergy(magBPA, trueMask);
artifactMF  = computeArtifactEnergy(magMF,  trueMask);
fprintf('ADDITIONAL metric (not from paper) - fraction of image energy OUTSIDE true target support:\n');
fprintf('  Classical BPA:            %.3f\n', artifactBPA);
fprintf('  Proposed measurement-fusion: %.3f\n', artifactMF);

% peak amplitude at each true scatterer location -> uniformity metric
peaksBPA = zeros(size(targets.pos,1), 1);
peaksMF  = zeros(size(targets.pos,1), 1);
for t = 1:size(targets.pos, 1)
    d = sqrt((X - targets.pos(t,1)).^2 + (Y - targets.pos(t,2)).^2);
    localMask = d < tol;
    if any(localMask(:))
        peaksBPA(t) = max(magBPA(localMask));
        peaksMF(t)  = max(magMF(localMask));
    end
end
uBPA = computeUniformity(peaksBPA);
uMF  = computeUniformity(peaksMF);
fprintf('ADDITIONAL metric (not from paper) - coefficient of variation of scatterer peak amplitudes:\n');
fprintf('  Classical BPA:              %.3f\n', uBPA.coefficientOfVariation);
fprintf('  Proposed measurement-fusion: %.3f\n', uMF.coefficientOfVariation);

save(fullfile(cfg.resultsDir, 'point_array_data.mat'), ...
    'targets', 'xAxis', 'yAxis', 'imgBPA', 'imgMF', 'artifactBPA', 'artifactMF', 'uBPA', 'uMF');

fprintf('Saved: point_array_BPA.png, point_array_proposed.png, point_array_comparison.png, point_array_data.mat\n');
end
