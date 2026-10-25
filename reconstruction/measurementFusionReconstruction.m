function [img, diag] = measurementFusionReconstruction(aperture, echo, k, voxelPos, cfg)
% MEASUREMENTFUSIONRECONSTRUCTION  Complete proposed algorithm, Eq. (13):
%
%   sigma_hat'(x,y,z) = w(x,y,z) * INT s(theta_uni,phi_uni,k) *
%                        exp(2j*k*RR) * delta_{alpha<Btheta/2} ...
%                        dtheta_uni dphi_uni dk
%
% which combines, in order:
%   1) Uniform sub-aperture construction               -- Eq. (4)-(7)
%   2) HPBW-gated back-projection over that sub-aperture -- Eq. (8),(10),(12)
%   3) Spatial-grid weighting (Hadamard product)         -- Eq. (9)-(11)
%
% -----------------------------------------------------------------------
% OUR EXTENSIONS (all OPTIONAL, all OFF by default -- see config.m and
% ADAPTIVE_METHODS.md). With the default configuration
%   cfg.weighting.mode = 'original'   and   cfg.gating.mode = 'fixed'
% this function is numerically identical to the original paper
% implementation.
%
%   cfg.weighting.mode : 'original' -> Eq.(9) geometry weight 1/count
%                        'adaptive' -> geometry weight modulated by a
%                                      bounded reconstructed-magnitude
%                                      confidence factor (computeAdaptiveWeight.m)
%   cfg.gating.mode    : 'fixed'    -> fixed nominal cfg.HPBW_deg gate
%                        'adaptive' -> per-sample effective HPBW from
%                                      sampling geometry (computeAdaptiveHPBW.m)
%
% The four supported combinations:
%   A  original weighting + fixed gating   (paper, DEFAULT)
%   B  adaptive weighting + fixed gating
%   C  original weighting + adaptive gating
%   D  adaptive weighting + adaptive gating
% -----------------------------------------------------------------------
%
%   aperture : struct from generateCurvedArray.m (RAW curved aperture)
%   echo     : Ns x Nf echo data over the RAW aperture (same sample
%              ordering as `aperture`)
%   k        : Nf x 1
%   voxelPos : Nv x 3
%   cfg      : configuration struct
%
%   img  : Nv x 1 complex reconstructed field, sigma_hat'
%   diag : struct with intermediate results. Original fields preserved:
%          subAp, weight, count, rawGatedImage. Added fields:
%          geometryWeight, signalConfidence, adaptiveFactor,
%          effectiveHPBW (summary), weightingMode, gatingMode,
%          runtimeSeconds, executionMode, chunkSizeUsed.

tStart = tic;

% --- defaults for a cfg built without the extension fields (legacy) ---
weightingMode = 'original';
if isfield(cfg, 'weighting') && isfield(cfg.weighting, 'mode')
    weightingMode = lower(cfg.weighting.mode);
end
gatingMode = 'fixed';
if isfield(cfg, 'gating') && isfield(cfg.gating, 'mode')
    gatingMode = lower(cfg.gating.mode);
end

% OUR EXTENSION: resolve performance mode so the chunkSize used here for
% the spatial weights matches the one used inside backProjectionCore.
cfgR = cfg;
if isfield(cfg, 'performance')
    % Ns for the sub-aperture is not yet known; use Nv as a proxy for the
    % auto chunk-size estimate here (backProjectionCore re-resolves it
    % exactly with the true Ns). chunk size only affects memory/speed,
    % never the result.
    cfgR = applyPerformanceMode(cfg, size(voxelPos,1), size(voxelPos,1), numel(k));
end
chunkSize = cfgR.recon.chunkSize;

subAp = createUniformSubAperture(aperture);           % Eq. (4)-(7)
echoUni = echo(subAp.idx, :);                          % Eq. (7) selection applied to the measured data

% ---- Eq. (12) gate: fixed nominal, or OUR adaptive per-sample HPBW ---
if strcmp(gatingMode, 'adaptive')
    [gateHPBW_deg, hpbwSummary] = computeAdaptiveHPBW(subAp, cfg);   % Ns x 1
else
    gateHPBW_deg = cfg.HPBW_deg;                                     % scalar (original)
    hpbwSummary = struct('mode', 'fixed', 'nominal_deg', cfg.HPBW_deg, ...
        'min_deg', cfg.HPBW_deg, 'max_deg', cfg.HPBW_deg, 'mean_deg', cfg.HPBW_deg);
end

% Eq. (9)-(10): geometry-only weight (echo/frequency independent),
% computed once. Uses the SAME per-sample HPBW as the gated
% back-projection so `count` is consistent.
geometryWeight = calculateSpatialWeights(subAp.pos, subAp.beamDir, voxelPos, gateHPBW_deg, chunkSize);

% Eq. (8) + Eq. (12): HPBW-gated back-projection over the uniform sub-aperture.
[rawImg, count, execInfo] = backProjectionCore(subAp.pos, subAp.beamDir, echoUni, k, voxelPos, cfg, true, gateHPBW_deg);

% ---- weighting: original Eq.(9) or OUR adaptive modulation ----------
if strcmp(weightingMode, 'adaptive')
    [weight, wdiag] = computeAdaptiveWeight(rawImg, count, geometryWeight, cfg);
    signalConfidence = wdiag.signalConfidence;
    adaptiveFactor   = wdiag.adaptiveFactor;
else
    weight = geometryWeight;                          % Eq. (9), original
    signalConfidence = [];
    adaptiveFactor   = [];
end

% Eq. (11): Hadamard product of the gated image with the weight.
img = rawImg .* weight;

% ---- diagnostics (original fields first, for backward compatibility) -
diag.subAp = subAp;
diag.weight = weight;
diag.count = count;
diag.rawGatedImage = rawImg;

diag.geometryWeight   = geometryWeight;
diag.signalConfidence = signalConfidence;
diag.adaptiveFactor   = adaptiveFactor;
diag.effectiveHPBW    = hpbwSummary;      % summary only (not Ns-length) to keep diag small
if strcmp(gatingMode, 'adaptive') && numel(gateHPBW_deg) <= 20000
    diag.effectiveHPBW.perSample_deg = gateHPBW_deg;
end
diag.weightingMode  = weightingMode;
diag.gatingMode     = gatingMode;
diag.executionMode  = execInfo.executionMode;
diag.chunkSizeUsed  = execInfo.chunkSize;
diag.runtimeSeconds = toc(tStart);

end
