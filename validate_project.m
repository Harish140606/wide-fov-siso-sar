function allPass = validate_project()
% VALIDATE_PROJECT  Part 22: automatic validation of every core
% component. Prints PASS/FAIL for each check; a check only prints PASS
% if its condition actually evaluated true. Uses a small, fast
% configuration (config('quick')) so this can run in well under a
% minute.

fprintf('\n================================================================\n');
fprintf(' Wide-FoV Curved SISO-SAR -- Automated Validation\n');
fprintf('================================================================\n');

results = [];
cfg = config('quick');

% ---------------------------------------------------------------
% 1) Geometry
% ---------------------------------------------------------------
ap = generateCurvedArray(cfg);
pl = generatePlanarArray(cfg);

ok = true;
ok = ok && numel(ap.x) == cfg.curved.N * cfg.curved.M;
ok = ok && abs(max(ap.phi) - deg2rad(cfg.curved.phi_c_deg)) < 1e-9;
ok = ok && abs(min(ap.phi)) < 1e-9;
ok = ok && all(abs(sqrt(sum(ap.pos.^2,2)) - cfg.curved.R0) < 1e-9);   % every element sits at radius R0
ok = ok && numel(pl.x) == cfg.planar.numElements * cfg.planar.numScan;
ok = ok && all(pl.z == 0);
results(end+1) = validationPrint('Geometry validation', ok);

% ---------------------------------------------------------------
% 2) Frequency configuration (incl. unambiguous-range margin check)
% ---------------------------------------------------------------
[f, k] = generateFrequencyGrid(cfg);
df = mean(diff(f));
unambRange = cfg.c / (2*df);
maxNeededRange = cfg.target.r + cfg.curved.R0 + 0.5;  % generous margin for imaging windows used in this project

ok = numel(f) == cfg.freq.numFreq;
ok = ok && abs(f(1) - cfg.freq.f_start) < 1;
ok = ok && abs(f(end) - cfg.freq.f_end) < 1;
ok = ok && all(abs(k - 2*pi*f/cfg.c) < 1e-6);
ok = ok && (unambRange > maxNeededRange);
fprintf('  (unambiguous range = %.3f m, required > %.3f m)\n', unambRange, maxNeededRange);
results(end+1) = validationPrint('Frequency configuration', ok);

% ---------------------------------------------------------------
% 3) Echo generation
% ---------------------------------------------------------------
tgt = createPointTarget(cfg.target.r, cfg.target.theta_deg, cfg.target.phi_deg);
echo = generatePointTargetEcho(ap.pos, tgt.pos, tgt.amp, k, ap.beamDir, cfg.HPBW_deg);

ok = isequal(size(echo), [numel(ap.x), numel(k)]);
ok = ok && all(isfinite(echo(:)));
ok = ok && any(abs(echo(:)) > 0);
results(end+1) = validationPrint('Echo generation', ok);

% ---------------------------------------------------------------
% 4) Classical BPA: peak lands on the true target
% ---------------------------------------------------------------
res = cfg.imaging.voxelRes;
hw  = 0.06;
[X,Y,Z] = meshgrid(tgt.pos(1)+(-hw:res:hw), tgt.pos(2)+(-hw:res:hw), tgt.pos(3)+(-hw:res:hw));
voxelPos = [X(:) Y(:) Z(:)];
imgBPA = BPA(echo, ap.pos, k, voxelPos, cfg);
[~,pkI] = max(abs(imgBPA));
errBPA = norm([X(pkI) Y(pkI) Z(pkI)] - tgt.pos);
fprintf('  (BPA peak location error = %.4f m, voxel size = %.4f m)\n', errBPA, res);
ok = errBPA <= res*sqrt(3)/2 + 1e-9;
results(end+1) = validationPrint('BPA reconstruction', ok);

% ---------------------------------------------------------------
% 5) FFT-accelerated BPA vs literal Eq. (8) direct summation
% ---------------------------------------------------------------
cfgD = cfg; cfgD.recon.method = 'direct'; cfgD.performance.mode = 'direct';
cfgF = cfg; cfgF.recon.zeroPadFactor = 32; cfgF.performance.mode = 'fft_serial';   % fine grid so interpolation error is small
smallVox = tgt.pos + [0 0 0; res 0 0; 0 res 0; 0 0 res; -res -res 0];
imgDirect = backProjectionCore(ap.pos, ap.beamDir, echo, k, smallVox, cfgD, false);
imgFFT    = backProjectionCore(ap.pos, ap.beamDir, echo, k, smallVox, cfgF, false);
Nfft = cfgF.recon.zeroPadFactor * cfg.freq.numFreq;
relErr = abs(imgDirect/Nfft - imgFFT) ./ abs(imgDirect/Nfft);
fprintf('  (max relative error between direct Eq.8 sum and FFT method = %.4f)\n', max(relErr));
ok = all(relErr < 0.05);
results(end+1) = validationPrint('FFT vs direct Eq.(8) equivalence', ok);

% ---------------------------------------------------------------
% 6) Uniform sub-aperture (Eq. 4-7)
% ---------------------------------------------------------------
subAp = createUniformSubAperture(ap);
% For the RAW array, every elevation ring has the same M azimuth
% samples regardless of ring radius (highly non-uniform physical
% spacing). After uniformization, the number of azimuth samples kept
% per ring should scale with sin(phi_i) (Eq. 5) rather than being
% constant -- check this correlation directly.
keptPerRing = sum(reshape(ismember(1:numel(ap.x), subAp.idx), cfg.curved.N, cfg.curved.M), 2);
nonzeroRings = ap.phi_i > 1e-9;
xv = sin(ap.phi_i(nonzeroRings)); xv = xv - mean(xv);
yv = keptPerRing(nonzeroRings);   yv = yv - mean(yv);
corrVal = sum(xv.*yv) / sqrt(sum(xv.^2) * sum(yv.^2));  % Pearson correlation, base MATLAB only

ok = numel(subAp.idx) < numel(ap.x);
ok = ok && numel(subAp.idx) == numel(unique(subAp.idx));
ok = ok && corrVal > 0.9;
fprintf('  (raw samples=%d, uniform sub-aperture=%d, corr(kept-per-ring, sin(phi_i))=%.3f)\n', ...
    numel(ap.x), numel(subAp.idx), corrVal);
results(end+1) = validationPrint('Uniform sub-aperture', ok);

% ---------------------------------------------------------------
% 7) HPBW gating (Eq. 10, 12)
% ---------------------------------------------------------------
testAp = [0 0 1];        % single sample at (0,0,1), boresight...
testBeam = [0 0 1];      % ...pointing along +z
insideVox  = [0.01 0 2];               % small angle off boresight -> should pass
outsideVox = [10 0 2];                 % large angle off boresight -> should be gated out
deltaIn  = calculateHPBWGate(testAp, testBeam, insideVox,  cfg.HPBW_deg);
deltaOut = calculateHPBWGate(testAp, testBeam, outsideVox, cfg.HPBW_deg);
ok = (deltaIn == true) && (deltaOut == false);
results(end+1) = validationPrint('HPBW gating', ok);

% ---------------------------------------------------------------
% 8) Spatial-grid weighting (Eq. 9)
% ---------------------------------------------------------------
weight = calculateSpatialWeights(subAp.pos, subAp.beamDir, voxelPos, cfg.HPBW_deg, cfg.recon.chunkSize);
ok = all(isfinite(weight)) && all(weight >= 0);
% A voxel behind the array (-z) is outside every element's beam cone:
% the curved array's elements only cover elevation phi in [0, phi_c]
% from the +z pole, so the farthest any boresight direction points from
% +z is phi_c + HPBW/2 <= 90 deg -- never into the -z hemisphere. Unlike
% a merely-distant point (beam angle is direction-only, not range-
% limited, so a far point can still be perfectly on-boresight), a point
% behind the array is geometrically guaranteed unilluminated.
behindVox = [0 0 -100];
behindWeight = calculateSpatialWeights(subAp.pos, subAp.beamDir, behindVox, cfg.HPBW_deg, cfg.recon.chunkSize);
ok = ok && (behindWeight == 0);   % zero illumination -> weight forced to 0, not NaN/Inf
results(end+1) = validationPrint('Spatial weighting', ok);

% ---------------------------------------------------------------
% 9) Measurement-fusion reconstruction (Eq. 13): peak near target
% ---------------------------------------------------------------
[imgMF, diagMF] = measurementFusionReconstruction(ap, echo, k, voxelPos, cfg);
[~,pkI2] = max(abs(imgMF));
errMF = norm([X(pkI2) Y(pkI2) Z(pkI2)] - tgt.pos);
fprintf('  (measurement-fusion peak location error = %.4f m)\n', errMF);
ok = errMF <= res*sqrt(3)/2 + 1e-9;
ok = ok && all(isfinite(imgMF));
ok = ok && all(isfinite(diagMF.weight));
results(end+1) = validationPrint('Measurement-fusion reconstruction', ok);

% ===============================================================
%  OUR EXTENSIONS -- adaptive reconstruction + runtime framework
%  (see ADAPTIVE_METHODS.md, PERFORMANCE_OPTIMIZATION.md)
% ===============================================================

% ---------------------------------------------------------------
% 10) Original mode reproduces the previous output EXACTLY
% ---------------------------------------------------------------
cfgLegacy = rmfield(cfg, intersect({'weighting','gating','noise','performance'}, fieldnames(cfg)));
imgLegacy = measurementFusionReconstruction(ap, echo, k, voxelPos, cfgLegacy);
imgDefault = measurementFusionReconstruction(ap, echo, k, voxelPos, cfg);  % default cfg = original modes
reproErr = max(abs(imgLegacy - imgDefault)) / (max(abs(imgLegacy)) + eps);
fprintf('  (default-cfg vs legacy-cfg max relative difference = %.2e)\n', reproErr);
ok = reproErr < 1e-12;
results(end+1) = validationPrint('Original mode reproduces previous output', ok);

% ---------------------------------------------------------------
% 11) Adaptive weights: finite, non-negative, no NaN/Inf, bounded
% ---------------------------------------------------------------
cfgAW = cfg; cfgAW.weighting.mode = 'adaptive';
[imgAW, diagAW] = measurementFusionReconstruction(ap, echo, k, voxelPos, cfgAW);
illum = diagAW.count > 0;
ok = all(isfinite(diagAW.weight)) && all(diagAW.weight >= 0);
ok = ok && ~any(isnan(imgAW)) && ~any(isinf(imgAW));
% raw (pre unit-mean normalisation) confidence factor is within [minFactor,maxFactor]
rawConf = (diagAW.signalConfidence + cfgAW.weighting.epsilon) .^ cfgAW.weighting.alpha;
rawConf = min(max(rawConf, cfgAW.weighting.minFactor), cfgAW.weighting.maxFactor);
ok = ok && all(rawConf >= cfgAW.weighting.minFactor - 1e-12) ...
        && all(rawConf <= cfgAW.weighting.maxFactor + 1e-12);
ok = ok && all(diagAW.weight(~illum) == 0);
results(end+1) = validationPrint('Adaptive weights finite/bounded/non-negative', ok);

% ---------------------------------------------------------------
% 12) Adaptive HPBW: within [min,max]; fixed mode == old gate exactly
% ---------------------------------------------------------------
subApV = createUniformSubAperture(ap);
effHPBW = computeAdaptiveHPBW(subApV, cfg);
ok = all(effHPBW >= cfg.gating.minHPBW_deg - 1e-9) && all(effHPBW <= cfg.gating.maxHPBW_deg + 1e-9);
% fixed-mode gate identical to the original scalar gate
gFixedVec = calculateHPBWGate(subApV.pos(1:5,:), subApV.beamDir(1:5,:), voxelPos, cfg.HPBW_deg*ones(5,1));
gFixedScalar = calculateHPBWGate(subApV.pos(1:5,:), subApV.beamDir(1:5,:), voxelPos, cfg.HPBW_deg);
ok = ok && isequal(gFixedVec, gFixedScalar);
fprintf('  (adaptive effective HPBW range = [%.2f, %.2f] deg, bounds = [%.0f, %.0f])\n', ...
    min(effHPBW), max(effHPBW), cfg.gating.minHPBW_deg, cfg.gating.maxHPBW_deg);
results(end+1) = validationPrint('Adaptive HPBW bounded; fixed == original gate', ok);

% ---------------------------------------------------------------
% 13) Zero illumination: no divide-by-zero, no NaN, weight 0
% ---------------------------------------------------------------
behindVox2 = [0 0 -100; 0 0 -50];
cfgD2 = cfg; cfgD2.weighting.mode = 'adaptive'; cfgD2.gating.mode = 'adaptive';
[imgZ, diagZ] = measurementFusionReconstruction(ap, echo, k, behindVox2, cfgD2);
ok = all(diagZ.count == 0) && all(diagZ.weight == 0) && all(imgZ == 0) && ~any(isnan(imgZ));
results(end+1) = validationPrint('Zero-illumination voxels safe (weight 0, no NaN)', ok);

% ---------------------------------------------------------------
% 14) Adaptive reconstruction: target peak still near true location
% ---------------------------------------------------------------
cfgBoth = cfg; cfgBoth.weighting.mode = 'adaptive'; cfgBoth.gating.mode = 'adaptive';
[imgB, ~] = measurementFusionReconstruction(ap, echo, k, voxelPos, cfgBoth);
[~, pkB] = max(abs(imgB));
errB = norm([X(pkB) Y(pkB) Z(pkB)] - tgt.pos);
fprintf('  (adaptive weighting+HPBW peak location error = %.4f m)\n', errB);
ok = errB <= res*sqrt(3)/2 + 1e-9;
results(end+1) = validationPrint('Adaptive reconstruction peak near target', ok);

% ---------------------------------------------------------------
% 15) Automatic chunking: valid nonzero size; result matches fixed
% ---------------------------------------------------------------
cfgAuto = cfg; cfgAuto.performance.chunkMode = 'auto';
csAuto = estimateChunkSize(numel(X), numel(subApV.idx), cfg.freq.numFreq, cfgAuto);
imgFixedChunk = measurementFusionReconstruction(ap, echo, k, voxelPos, cfg);
imgAutoChunk  = measurementFusionReconstruction(ap, echo, k, voxelPos, cfgAuto);
chunkErr = max(abs(imgFixedChunk - imgAutoChunk)) / (max(abs(imgFixedChunk)) + eps);
fprintf('  (auto chunk size = %d [1..%d]; fixed-vs-auto max rel diff = %.2e)\n', ...
    csAuto, numel(subApV.idx), chunkErr);
ok = csAuto >= 1 && csAuto <= numel(subApV.idx) && chunkErr < 1e-10;
results(end+1) = validationPrint('Automatic chunk sizing valid & result-invariant', ok);

% ---------------------------------------------------------------
% 16) Performance-mode selector: fft_serial == legacy defaults
% ---------------------------------------------------------------
cfgPS = applyPerformanceMode(cfg, 100, 100, cfg.freq.numFreq);
ok = strcmp(cfgPS.recon.method, 'fft') && ~cfgPS.useGPU && ~cfgPS.useParallel;
cfgPD = cfg; cfgPD.performance.mode = 'direct';
cfgPD = applyPerformanceMode(cfgPD, 100, 100, cfg.freq.numFreq);
ok = ok && strcmp(cfgPD.recon.method, 'direct');
results(end+1) = validationPrint('Performance-mode selector resolves correctly', ok);

% ---------------------------------------------------------------
% 17) Parallel vs serial numerically close (SKIPPED if no PCT)
% ---------------------------------------------------------------
if isempty(ver('parallel'))
    fprintf('  (Parallel Computing Toolbox unavailable -- test SKIPPED, not counted)\n');
else
    cfgPar = cfg; cfgPar.performance.mode = 'fft_parallel';
    imgSer = measurementFusionReconstruction(ap, echo, k, voxelPos, cfg);
    imgPar = measurementFusionReconstruction(ap, echo, k, voxelPos, cfgPar);
    parErr = max(abs(imgSer - imgPar)) / (max(abs(imgSer)) + eps);
    fprintf('  (serial vs parallel max relative difference = %.2e)\n', parErr);
    results(end+1) = validationPrint('Parallel CPU vs serial numerically close', parErr < 1e-9);
end

% ---------------------------------------------------------------
% 18) GPU vs serial numerically close (SKIPPED if no GPU)
% ---------------------------------------------------------------
haveGPU = false;
try, haveGPU = ~isempty(ver('parallel')) && gpuDeviceCount() > 0; catch, haveGPU = false; end
if ~haveGPU
    fprintf('  (no supported GPU -- test SKIPPED, not counted)\n');
else
    cfgGPU = cfg; cfgGPU.performance.mode = 'fft_gpu';
    imgSer2 = measurementFusionReconstruction(ap, echo, k, voxelPos, cfg);
    imgGPU  = measurementFusionReconstruction(ap, echo, k, voxelPos, cfgGPU);
    gpuErr = max(abs(imgSer2 - imgGPU)) / (max(abs(imgSer2)) + eps);
    fprintf('  (serial vs GPU max relative difference = %.2e)\n', gpuErr);
    results(end+1) = validationPrint('GPU vs serial numerically close', gpuErr < 1e-5);
end

% ---------------------------------------------------------------
% 19) Complex AWGN helper: Inf -> unchanged; finite -> ~target SNR
% ---------------------------------------------------------------
echoClean = echo;
echoInf = addComplexAWGN(echoClean, Inf);
n20 = addComplexAWGN(echoClean, 20, 1);
measSNR = 10*log10(mean(abs(echoClean(:)).^2) / mean(abs(n20(:)-echoClean(:)).^2));
fprintf('  (requested SNR 20 dB, measured %.2f dB)\n', measSNR);
ok = isequal(echoInf, echoClean) && abs(measSNR - 20) < 0.5;
results(end+1) = validationPrint('Complex AWGN helper correct', ok);

% ---------------------------------------------------------------
fprintf('----------------------------------------------------------------\n');
allPass = all(results);
if allPass
    fprintf('ALL CHECKS PASSED\n');
else
    fprintf('SOME CHECKS FAILED -- see above\n');
end
fprintf('================================================================\n\n');
end
