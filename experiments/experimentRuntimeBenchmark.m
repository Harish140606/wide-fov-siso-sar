function out = experimentRuntimeBenchmark(cfg)
% EXPERIMENTRUNTIMEBENCHMARK  OUR EXTENSION (implementation contribution).
%
% Objective, reproducible benchmark of the reconstruction execution
% strategies. NOTE: FFT evaluation, chunking, parfor and gpuArray
% support ALREADY EXISTED in backProjectionCore.m -- this experiment is
% the quantitative comparison / auto-selection layer, not a new radar
% algorithm. See PERFORMANCE_OPTIMIZATION.md.
%
% Methods (unavailable ones are marked 'unavailable', never faked):
%   1. Direct BPA                     (literal Eq.8 sum)
%   2. FFT serial
%   3. FFT serial + auto chunking
%   4. FFT parallel CPU               (needs Parallel Computing Toolbox)
%   5. FFT GPU                        (needs supported GPU)
%
% Accuracy reference: the DIRECT method on the SAME small problem
% (images peak-normalised before comparison so the known Nfft scale
% factor between 'direct' and 'fft' does not appear as an error).
%
% For each method: runtimeSeconds, speedupVsDirect, speedupVsFFTSerial,
% peakLocationError, relativeImageErrorVsReference, memoryEstimateMB,
% actualExecutionMode, chunkSizeUsed.

if nargin < 1 || isempty(cfg)
    cfg = config('quick');
end

fprintf('\n==================================================================\n');
fprintf(' Runtime Benchmark  (mode: %s)\n', cfg.mode);
fprintf('==================================================================\n');

% ---- small, direct-feasible problem: single point target, 3-D window
ap = generateCurvedArray(cfg);
[~, k] = generateFrequencyGrid(cfg);
tgt = createPointTarget(cfg.target.r, cfg.target.theta_deg, cfg.target.phi_deg);
echo = generatePointTargetEcho(ap.pos, tgt.pos, tgt.amp, k, ap.beamDir, cfg.HPBW_deg);
echo = maybeAddNoise(echo, cfg);

res = cfg.imaging.voxelRes;
hw  = 0.06;
[X,Y,Z] = meshgrid(tgt.pos(1)+(-hw:res:hw), tgt.pos(2)+(-hw:res:hw), tgt.pos(3)+(-hw:res:hw));
voxelPos = [X(:) Y(:) Z(:)];
Nv = size(voxelPos,1); Ns = size(ap.pos,1); Nf = numel(k);
fprintf('Problem: Ns=%d aperture samples, Nf=%d freqs, Nv=%d voxels\n', Ns, Nf, Nv);

hasPCT = ~isempty(ver('parallel'));
haveGPU = false;
try, haveGPU = hasPCT && gpuDeviceCount() > 0; catch, haveGPU = false; end

specs = {
    'Direct BPA',                 'direct',       'fixed'
    'FFT serial',                 'fft_serial',   'fixed'
    'FFT serial + auto chunking', 'fft_serial',   'auto'
    'FFT parallel CPU',           'fft_parallel', 'fixed'
    'FFT GPU',                    'fft_gpu',      'fixed'
};
nM = size(specs,1);

R = struct();
R.method               = specs(:,1)';
R.available            = true(1,nM);
R.runtimeSeconds       = nan(1,nM);
R.speedupVsDirect      = nan(1,nM);
R.speedupVsFFTSerial   = nan(1,nM);
R.peakLocationError    = nan(1,nM);
R.relativeImageError   = nan(1,nM);
R.memoryEstimateMB     = nan(1,nM);
R.actualExecutionMode  = repmat({'unavailable'},1,nM);
R.chunkSizeUsed        = nan(1,nM);

R.available(4) = hasPCT;
R.available(5) = haveGPU;

% reference image (direct)
cfgRef = cfg; cfgRef.performance.mode = 'direct';
imgRef = backProjectionCore(ap.pos, ap.beamDir, echo, k, voxelPos, cfgRef, false);
imgRefN = imgRef / (max(abs(imgRef)) + eps);

for m = 1:nM
    if ~R.available(m)
        fprintf('  %-28s : unavailable (skipped, not benchmarked)\n', specs{m,1});
        continue;
    end
    c = cfg;
    c.performance.mode = specs{m,2};
    c.performance.chunkMode = specs{m,3};

    [~, ~, info0] = backProjectionCore(ap.pos, ap.beamDir, echo, k, voxelPos, c, false); % warm-up / resolve
    t = tic;
    [img, ~, info] = backProjectionCore(ap.pos, ap.beamDir, echo, k, voxelPos, c, false);
    R.runtimeSeconds(m) = toc(t);

    imgN = img / (max(abs(img)) + eps);
    R.relativeImageError(m) = norm(imgN - imgRefN) / (norm(imgRefN) + eps);
    [~, pk] = max(abs(img));
    R.peakLocationError(m) = norm([X(pk) Y(pk) Z(pk)] - tgt.pos);
    R.actualExecutionMode{m} = info.executionMode;
    R.chunkSizeUsed(m) = info.chunkSize;
    [~, estMB] = estimateChunkSize(Nv, Ns, Nf, c);
    R.memoryEstimateMB(m) = estMB;

    fprintf('  %-28s : %7.3f s | mode=%s | chunk=%d | relErr=%.2e | pkErr=%.4f m\n', ...
        specs{m,1}, R.runtimeSeconds(m), R.actualExecutionMode{m}, ...
        R.chunkSizeUsed(m), R.relativeImageError(m), R.peakLocationError(m));
end

tDirect = R.runtimeSeconds(1);
tFFT    = R.runtimeSeconds(2);
R.speedupVsDirect    = tDirect ./ R.runtimeSeconds;
R.speedupVsFFTSerial = tFFT    ./ R.runtimeSeconds;
R.speedupVsDirect(~R.available)    = NaN;
R.speedupVsFFTSerial(~R.available) = NaN;

% ---- save .mat / .csv ----
save(fullfile(cfg.resultsDir, 'runtime_benchmark.mat'), 'R');
writeBenchmarkCSV(fullfile(cfg.resultsDir, 'runtime_benchmark.csv'), R);

% ---- plots (only available methods) ----
av = R.available;
names = R.method(av);

barPlot(names, R.runtimeSeconds(av), 'runtime / s', ...
    'Runtime vs Method', fullfile(cfg.resultsDir, 'runtime_comparison.png'));
barPlot(names, R.speedupVsDirect(av), 'speedup vs Direct BPA (x)', ...
    'Speedup vs Method', fullfile(cfg.resultsDir, 'speedup_comparison.png'));
barPlot(names, R.relativeImageError(av), 'relative image error vs direct reference', ...
    'Reconstruction Error vs Reference (peak-normalised)', ...
    fullfile(cfg.resultsDir, 'reconstruction_error_comparison.png'));

out = R;
fprintf('\nSaved: runtime_benchmark.mat, runtime_benchmark.csv,\n');
fprintf('       runtime_comparison.png, speedup_comparison.png, reconstruction_error_comparison.png\n');
end

% ------------------------------------------------------------------
function writeBenchmarkCSV(path, R)
fid = fopen(path, 'w');
fprintf(fid, ['method,available,runtimeSeconds,speedupVsDirect,speedupVsFFTSerial,', ...
              'peakLocationError,relativeImageErrorVsReference,memoryEstimateMB,', ...
              'actualExecutionMode,chunkSizeUsed\n']);
for m = 1:numel(R.method)
    if R.available(m), avStr = 'yes'; else, avStr = 'unavailable'; end
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
        R.method{m}, avStr, ...
        num2str(R.runtimeSeconds(m)), num2str(R.speedupVsDirect(m)), ...
        num2str(R.speedupVsFFTSerial(m)), num2str(R.peakLocationError(m)), ...
        num2str(R.relativeImageError(m)), num2str(R.memoryEstimateMB(m)), ...
        R.actualExecutionMode{m}, num2str(R.chunkSizeUsed(m)));
end
fclose(fid);
end

% ------------------------------------------------------------------
function barPlot(names, vals, ylab, ttl, savePath)
fig = figure('Visible','on','Position',[60 60 720 420]);
drawnow;
bar(vals);
set(gca, 'XTickLabel', names, 'XTickLabelRotation', 20, 'FontSize', 9);
ylabel(ylab); title(ttl, 'Interpreter', 'none'); grid on;
if ~isempty(savePath)
    exportgraphics(fig, savePath, 'Resolution', 150);
end
end
