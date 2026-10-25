function out = experimentNoiseRobustness(cfg)
% EXPERIMENTNOISEROBUSTNESS  OUR EXTENSION. Because the adaptive weighting
% uses reconstructed signal magnitude, it could in principle amplify
% noise. This experiment quantifies that directly by sweeping SNR and
% comparing:
%
%   - Original Measurement Fusion            (fixed HPBW + 1/count)
%   - Adaptive Weighting                     (fixed HPBW + adaptive weight)
%   - Adaptive Weighting + Adaptive HPBW
%
% SNR sweep: [Inf, 30, 20, 10] dB  (complex AWGN, utils/addComplexAWGN.m)
% Scene: planar point array (known true support -> clean metrics).
%
% Metrics (all ADDITIONAL - not reported in original paper):
%   - artifact energy fraction (energy outside true support)
%   - uniformity: coefficient of variation of true-scatterer peaks
%   - mean peak localization error (m)
%   - reconstruction quality: contrast (dB) = true-support mean vs background mean
%
% If an adaptive method degrades faster with noise than the original,
% the numbers are printed and saved unchanged.

if nargin < 1 || isempty(cfg)
    cfg = config('quick');
end

fprintf('\n==================================================================\n');
fprintf(' Noise / SNR Robustness  (mode: %s)\n', cfg.mode);
fprintf('==================================================================\n');

SNRlist = [Inf, 30, 20, 10];
methodNames = {'Original Measurement Fusion', 'Adaptive Weighting', 'Adaptive Weighting + Adaptive HPBW'};
variants = { struct('w','original','g','fixed'), ...
             struct('w','adaptive','g','fixed'), ...
             struct('w','adaptive','g','adaptive') };

ap = generateCurvedArray(cfg);
[~, k] = generateFrequencyGrid(cfg);
z0 = 1.0; halfWidth = 0.4; spacing = 0.1;
targets = createPointArray(z0, halfWidth, spacing);
echoClean = generateSceneEcho(ap.pos, targets, k, ap.beamDir, cfg.HPBW_deg);

res = cfg.imaging.voxelRes; hw = 0.5;
xAxis = -hw:res:hw; yAxis = -hw:res:hw;
[X, Y] = meshgrid(xAxis, yAxis);
voxelPos = [X(:), Y(:), z0*ones(numel(X),1)];

tol = res * 1.5;
trueMask = false(size(X));
for t = 1:size(targets.pos,1)
    d = sqrt((X - targets.pos(t,1)).^2 + (Y - targets.pos(t,2)).^2);
    trueMask = trueMask | (d < tol);
end

nS = numel(SNRlist); nM = numel(methodNames);
artifact = nan(nS, nM);
peakCoV  = nan(nS, nM);
locErr   = nan(nS, nM);
contrast = nan(nS, nM);

for si = 1:nS
    echo = addComplexAWGN(echoClean, SNRlist(si), cfg.noise.seed);
    for mi = 1:nM
        c = cfg;
        c.weighting.mode = variants{mi}.w;
        c.gating.mode    = variants{mi}.g;
        img = measurementFusionReconstruction(ap, echo, k, voxelPos, c);
        mag = reshape(abs(img), size(X));

        artifact(si,mi) = computeArtifactEnergy(mag, trueMask);
        peaks = zeros(size(targets.pos,1),1);
        le = zeros(size(targets.pos,1),1);
        for t = 1:size(targets.pos,1)
            d = sqrt((X - targets.pos(t,1)).^2 + (Y - targets.pos(t,2)).^2);
            lm = d < 2*tol;
            sub = mag; sub(~lm) = -Inf;
            [pk, iMax] = max(sub(:));
            peaks(t) = pk;
            le(t) = sqrt((X(iMax)-targets.pos(t,1)).^2 + (Y(iMax)-targets.pos(t,2)).^2);
        end
        u = computeUniformity(peaks);
        peakCoV(si,mi) = u.coefficientOfVariation;
        locErr(si,mi)  = mean(le);
        contrast(si,mi) = 20*log10((mean(mag(trueMask))+eps)/(mean(mag(~trueMask))+eps));
    end
end

fprintf('\n%-8s %-36s %10s %8s %9s %9s\n', 'SNR_dB', 'Method', 'artifact', 'peakCoV', 'locErr_m', 'contr_dB');
for si = 1:nS
    for mi = 1:nM
        fprintf('%-8s %-36s %10.3f %8.3f %9.4f %9.2f\n', num2str(SNRlist(si)), methodNames{mi}, ...
            artifact(si,mi), peakCoV(si,mi), locErr(si,mi), contrast(si,mi));
    end
end
fprintf(['\n  All metrics above: Additional metric - not reported in original paper\n' ...
         '  Interpretation: if "Adaptive Weighting" artifact/CoV grows faster than\n' ...
         '  "Original Measurement Fusion" as SNR drops, adaptive weighting is\n' ...
         '  amplifying noise -- reported honestly either way.\n']);

curveFigure(SNRlist, artifact, methodNames, 'artifact energy fraction', ...
    'Noise Robustness - Artifact Energy vs SNR', fullfile(cfg.resultsDir, 'noise_artifact_vs_snr.png'));
curveFigure(SNRlist, peakCoV, methodNames, 'peak amplitude CoV (uniformity)', ...
    'Noise Robustness - Uniformity vs SNR', fullfile(cfg.resultsDir, 'noise_uniformity_vs_snr.png'));
curveFigure(SNRlist, locErr, methodNames, 'mean peak localization error / m', ...
    'Noise Robustness - Localization Error vs SNR', fullfile(cfg.resultsDir, 'noise_localization_vs_snr.png'));

out = struct('SNRlist', SNRlist, 'methodNames', {methodNames}, ...
             'artifact', artifact, 'peakCoV', peakCoV, 'locErr', locErr, 'contrast', contrast);
save(fullfile(cfg.resultsDir, 'noise_robustness_data.mat'), 'out');
fprintf('\nSaved: noise_artifact_vs_snr.png, noise_uniformity_vs_snr.png,\n');
fprintf('       noise_localization_vs_snr.png, noise_robustness_data.mat\n');
end

% ------------------------------------------------------------------
function curveFigure(SNRlist, M, names, ylab, ttl, savePath)
% x axis = ordinal index (SNRlist order: Inf,30,20,10), labelled with SNR
xi = 1:numel(SNRlist);
labs = cell(1,numel(SNRlist));
for j = 1:numel(SNRlist)
    if isinf(SNRlist(j)), labs{j} = 'Inf (clean)'; else, labs{j} = sprintf('%g', SNRlist(j)); end
end
fig = figure('Visible','on'); drawnow;
hold on; grid on; box on;
markers = {'-o','-s','-^','-d'};
for i = 1:size(M,2)
    plot(xi, M(:,i), markers{mod(i-1,4)+1}, 'LineWidth', 1.4, 'MarkerSize', 6);
end
set(gca, 'XTick', xi, 'XTickLabel', labs, 'XLim', [0.7 numel(xi)+0.3]);
xlabel('echo SNR / dB'); ylabel(ylab); title(ttl, 'Interpreter', 'none');
legend(names, 'Location', 'best', 'Interpreter', 'none');
if ~isempty(savePath), exportgraphics(fig, savePath, 'Resolution', 150); end
end
