function out = experimentAdaptiveMethods(cfg)
% EXPERIMENTADAPTIVEMETHODS  OUR EXTENSION. Head-to-head comparison of
% the original paper reconstruction and our adaptive variants:
%
%   1. Classical BPA                          (raw aperture, un-gated)
%   2. Original Measurement Fusion            (Eq.13: fixed HPBW + 1/count)
%   3. Adaptive Weighting                     (fixed HPBW + adaptive weight)
%   4. Adaptive HPBW                          (adaptive HPBW + 1/count)
%   5. Adaptive Weighting + Adaptive HPBW     (both)
%
% Scenes:
%   A. Planar point-array (Fig. 6 scene) -- full metric table
%   B. Resolution vs elevation angle      -- -3 dB azimuth width per method
%   C. Human point cloud                  -- only if cfg.mode == 'quick'
%                                            is NOT set to skip (runtime)
%
% Every metric that is NOT reported in the paper is printed with the tag
% "Additional metric - not reported in original paper".
%
% If an adaptive method does WORSE than the original on a metric, the
% number is printed and saved unchanged -- nothing is hidden.

if nargin < 1 || isempty(cfg)
    cfg = config('quick');
end

fprintf('\n==================================================================\n');
fprintf(' Adaptive Methods Comparison  (mode: %s)\n', cfg.mode);
fprintf('==================================================================\n');

methodNames = { ...
    'Classical BPA', ...
    'Original Measurement Fusion', ...
    'Adaptive Weighting', ...
    'Adaptive HPBW', ...
    'Adaptive Weighting + Adaptive HPBW'};

% cfg variants (weighting.mode / gating.mode) for methods 2..5
variants = { ...
    struct('weighting','original','gating','fixed'), ...
    struct('weighting','adaptive','gating','fixed'), ...
    struct('weighting','original','gating','adaptive'), ...
    struct('weighting','adaptive','gating','adaptive')};

ap = generateCurvedArray(cfg);
[~, k] = generateFrequencyGrid(cfg);

% =================================================================
% SCENE A: planar point array
% =================================================================
fprintf('\n--- Scene A: planar point array (z = 1 m) ---\n');
z0 = 1.0; halfWidth = 0.4; spacing = 0.1;
targets = createPointArray(z0, halfWidth, spacing);

echo = generateSceneEcho(ap.pos, targets, k, ap.beamDir, cfg.HPBW_deg);
echo = maybeAddNoise(echo, cfg);

hw = 0.5; res = cfg.imaging.voxelRes;
xAxis = -hw:res:hw; yAxis = -hw:res:hw;
[X, Y] = meshgrid(xAxis, yAxis);
voxelPos = [X(:), Y(:), z0*ones(numel(X),1)];

imgs = cell(1, 5);
imgs{1} = BPA(echo, ap.pos, k, voxelPos, cfg);
for m = 2:5
    c = cfg;
    c.weighting.mode = variants{m-1}.weighting;
    c.gating.mode    = variants{m-1}.gating;
    imgs{m} = measurementFusionReconstruction(ap, echo, k, voxelPos, c);
end

% --- metrics ---
tol = res * 1.5;
trueMask = false(size(X));
for t = 1:size(targets.pos,1)
    d = sqrt((X - targets.pos(t,1)).^2 + (Y - targets.pos(t,2)).^2);
    trueMask = trueMask | (d < tol);
end

A = struct();
A.method = methodNames;
A.artifactEnergyFraction = zeros(1,5);
A.peakCoV = zeros(1,5);
A.peakLocErrorMean_m = zeros(1,5);
A.contrast_dB = zeros(1,5);
A.peakSidelobe_dB = zeros(1,5);

for m = 1:5
    mag = reshape(abs(imgs{m}), size(X));
    A.artifactEnergyFraction(m) = computeArtifactEnergy(mag, trueMask);

    peaks = zeros(size(targets.pos,1),1);
    locErr = zeros(size(targets.pos,1),1);
    for t = 1:size(targets.pos,1)
        d = sqrt((X - targets.pos(t,1)).^2 + (Y - targets.pos(t,2)).^2);
        localMask = d < 2*tol;
        sub = mag; sub(~localMask) = -Inf;
        [pk, iMax] = max(sub(:));
        peaks(t) = pk;
        locErr(t) = sqrt((X(iMax)-targets.pos(t,1)).^2 + (Y(iMax)-targets.pos(t,2)).^2);
    end
    u = computeUniformity(peaks);
    A.peakCoV(m) = u.coefficientOfVariation;
    A.peakLocErrorMean_m(m) = mean(locErr);

    % contrast: mean true-support magnitude vs mean background magnitude (dB)
    A.contrast_dB(m) = 20*log10( (mean(mag(trueMask)) + eps) / (mean(mag(~trueMask)) + eps) );

    % peak sidelobe: strongest background pixel outside a guard band, dB below peak
    guard = false(size(X));
    for t = 1:size(targets.pos,1)
        d = sqrt((X - targets.pos(t,1)).^2 + (Y - targets.pos(t,2)).^2);
        guard = guard | (d < 3*tol);
    end
    A.peakSidelobe_dB(m) = 20*log10( (max(mag(~guard)) + eps) / (max(mag(:)) + eps) );
end

fprintf('\n%-38s %10s %8s %9s %9s %9s\n', 'Method', 'artifact', 'peakCoV', 'locErr_m', 'contr_dB', 'PSL_dB');
for m = 1:5
    fprintf('%-38s %10.3f %8.3f %9.4f %9.2f %9.2f\n', methodNames{m}, ...
        A.artifactEnergyFraction(m), A.peakCoV(m), A.peakLocErrorMean_m(m), ...
        A.contrast_dB(m), A.peakSidelobe_dB(m));
end
fprintf(['\n  artifact energy fraction  : Additional metric - not reported in original paper\n' ...
         '  peak amplitude CoV        : Additional metric - not reported in original paper\n' ...
         '  peak localization error   : Additional metric - not reported in original paper\n' ...
         '  contrast (dB)             : Additional metric - not reported in original paper\n' ...
         '  peak sidelobe level (dB)  : Additional metric - not reported in original paper\n']);

% honest verdict lines vs "Original Measurement Fusion" (m=2)
fprintf('\n  Change vs Original Measurement Fusion (negative artifact/CoV = better):\n');
for m = 3:5
    fprintf('   %-36s  d(artifact)=%+.3f  d(CoV)=%+.3f\n', methodNames{m}, ...
        A.artifactEnergyFraction(m) - A.artifactEnergyFraction(2), ...
        A.peakCoV(m) - A.peakCoV(2));
end

% --- figures ---
dR = cfg.imaging.dB_range;
imgCellDB = cellfun(@(v) reshape(toDB(v, dR), size(X)), imgs, 'UniformOutput', false);
plotComparison(imgCellDB, methodNames, xAxis, yAxis, 'x / m', 'y / m', dR, ...
    'Adaptive Methods - Planar Point Array', ...
    fullfile(cfg.resultsDir, 'adaptive_point_array_comparison.png'));

barMetricFigure(methodNames, ...
    {A.artifactEnergyFraction, A.peakCoV, A.peakLocErrorMean_m}, ...
    {'artifact energy fraction', 'peak amplitude CoV', 'mean peak loc error / m'}, ...
    'Adaptive Methods - Point-Array Metrics (all ADDITIONAL, not in paper)', ...
    fullfile(cfg.resultsDir, 'adaptive_point_array_metrics.png'));

% =================================================================
% SCENE B: resolution vs elevation angle (azimuth -3 dB width)
% =================================================================
fprintf('\n--- Scene B: azimuth resolution vs elevation angle ---\n');
phiList = cfg.resSweep.phi_deg_list;
r = cfg.resSweep.r; theta_deg = cfg.resSweep.theta_deg;
lineHW = 0.06; lineRes = 0.001;
azRes = nan(numel(phiList), 5);
for a = 1:numel(phiList)
    tgt = createPointTarget(r, theta_deg, phiList(a));
    th = deg2rad(theta_deg);
    azHat = [-sin(th), cos(th), 0];
    tLine = (-lineHW:lineRes:lineHW)';
    azLine = tgt.pos + tLine*azHat;
    echoC = generatePointTargetEcho(ap.pos, tgt.pos, tgt.amp, k, ap.beamDir, cfg.HPBW_deg);
    echoC = maybeAddNoise(echoC, cfg);
    azRes(a,1) = computeResolution(tLine, abs(BPA(echoC, ap.pos, k, azLine, cfg)));
    for m = 2:5
        c = cfg; c.weighting.mode = variants{m-1}.weighting; c.gating.mode = variants{m-1}.gating;
        azRes(a,m) = computeResolution(tLine, abs(measurementFusionReconstruction(ap, echoC, k, azLine, c)));
    end
end
plotResolution(phiList, num2cell(azRes,1), methodNames, ...
    'Adaptive Methods - Azimuth Resolution vs Elevation Angle', 'resolution / cm', ...
    fullfile(cfg.resultsDir, 'adaptive_resolution_azimuth.png'));

B = struct('phiList', phiList, 'azRes_m', azRes);
B.method = methodNames;

% =================================================================
% SCENE C: human point cloud (only when runtime is practical)
% =================================================================
C = [];
runHuman = ~(isfield(cfg,'adaptive') && isfield(cfg.adaptive,'skipHuman') && cfg.adaptive.skipHuman);
if runHuman
    fprintf('\n--- Scene C: human point cloud (y = 0 m slice) ---\n');
    try
        human = createHumanPointCloud(0, 1.0, 1);
        xAxisH = -1.0:res:0.95;
        zAxisH = 1.0 + (-0.3:res:0.3);
        [XH, ZH] = meshgrid(xAxisH, zAxisH);
        voxelH = [XH(:), zeros(numel(XH),1), ZH(:)];
        echoH = generateHumanPointCloudEcho(ap.pos, human, k, ap.beamDir, cfg.HPBW_deg);
        echoH = maybeAddNoise(echoH, cfg);
        imgsH = cell(1,5);
        imgsH{1} = BPA(echoH, ap.pos, k, voxelH, cfg);
        for m = 2:5
            c = cfg; c.weighting.mode = variants{m-1}.weighting; c.gating.mode = variants{m-1}.gating;
            imgsH{m} = measurementFusionReconstruction(ap, echoH, k, voxelH, c);
        end
        imgCellH = cellfun(@(v) reshape(toDB(v, dR), size(XH)), imgsH, 'UniformOutput', false);
        plotComparison(imgCellH, methodNames, xAxisH, zAxisH, 'x / m (height)', 'z / m (depth)', dR, ...
            'Adaptive Methods - Human Point Cloud (y = 0 m)', ...
            fullfile(cfg.resultsDir, 'adaptive_human_comparison.png'));
        C = struct('xAxis', xAxisH, 'zAxis', zAxisH);
        C.method = methodNames;
        for m = 1:5, C.(sprintf('img%d', m)) = imgsH{m}; end
    catch ME
        fprintf('  Human scene skipped (%s)\n', ME.message);
    end
end

out.methodNames = methodNames;
out.sceneA = A;
out.sceneB = B;
out.sceneC = C;
save(fullfile(cfg.resultsDir, 'adaptive_methods_data.mat'), 'out');
fprintf('\nSaved: adaptive_point_array_comparison.png, adaptive_point_array_metrics.png,\n');
fprintf('       adaptive_resolution_azimuth.png, adaptive_methods_data.mat\n');
end

% ------------------------------------------------------------------
function barMetricFigure(methodNames, metricVals, metricLabels, ttl, savePath)
fig = figure('Visible','on','Position',[50 50 400*numel(metricVals) 380]);
drawnow;
tl = tiledlayout(fig, 1, numel(metricVals), 'Padding','compact','TileSpacing','compact');
for i = 1:numel(metricVals)
    nexttile(tl);
    bar(metricVals{i});
    set(gca, 'XTickLabel', methodNames, 'XTickLabelRotation', 30, 'FontSize', 8);
    ylabel(metricLabels{i});
    grid on;
end
title(tl, ttl, 'Interpreter', 'none');
if nargin >= 5 && ~isempty(savePath)
    exportgraphics(fig, savePath, 'Resolution', 150);
end
end
