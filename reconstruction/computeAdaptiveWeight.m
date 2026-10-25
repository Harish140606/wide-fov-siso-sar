function [weight, d] = computeAdaptiveWeight(rawImg, count, geometryWeight, cfg)
% COMPUTEADAPTIVEWEIGHT  OUR PROPOSED EXTENSION (NOT from the paper).
%
% The paper's Eq.(9) weight is geometry only:  w = 1 / count  (count =
% number of sub-aperture samples that illuminate the voxel). It cannot
% tell a genuine, consistently-reconstructed scatterer from a
% low-amplitude interference / clutter voxel that happens to have the
% same illumination count.
%
% We keep that geometry weight and modulate it by a BOUNDED confidence
% factor derived from the already-computed gated back-projection
% magnitude |rawImg| (no extra reconstruction pass):
%
%   geometryWeight = 1 ./ max(count,1)                  (== Eq.9)
%   signalNorm     = |rawImg| ./ (max(|rawImg|) + eps)  in [0, ~1]
%   confidence     = (signalNorm + eps) .^ alpha
%   adaptiveFactor = clip(confidence, minFactor, maxFactor)
%   adaptiveFactor = adaptiveFactor ./ mean(adaptiveFactor | illuminated)
%   weight         = geometryWeight .* adaptiveFactor
%   weight(count==0) = 0
%
% Stability / honesty notes:
%   * epsilon floors the normalised magnitude so a near-zero voxel gives
%     confidence ~ eps^alpha, not 0 or Inf.
%   * clip(minFactor,maxFactor) prevents BOTH runaway boosting of weak
%     noisy voxels (upper clip) and total suppression of real but dim
%     scatterers (lower clip). Keep maxFactor modest (~2).
%   * the unit-mean normalisation over illuminated voxels keeps the
%     overall image scale comparable to the original weighting, so dB
%     displays and energy-fraction metrics stay directly comparable.
%   * alpha = 0  ->  confidence == 1 everywhere  ->  weight == Eq.9.
%   * This formula is heuristic and is NOT claimed to come from the
%     paper. See ADAPTIVE_METHODS.md and the SNR robustness experiment.
%
%   rawImg         : Nv x 1 complex, HPBW-gated back-projection
%   count          : Nv x 1, illumination counts
%   geometryWeight : Nv x 1, the Eq.(9) weight (1./count, 0 where count==0)
%   cfg            : configuration struct (cfg.weighting.*)
%
%   weight : Nv x 1, finite, non-negative, 0 where count==0
%   d      : diagnostics struct (signalConfidence, adaptiveFactor)

alpha     = cfg.weighting.alpha;
minFactor = cfg.weighting.minFactor;
maxFactor = cfg.weighting.maxFactor;
epsW      = cfg.weighting.epsilon;

mag = abs(rawImg(:));
mmax = max(mag);
if mmax < eps
    mmax = eps;
end

signalNorm = mag ./ (mmax + epsW);                  % in [0, ~1]
confidence = (signalNorm + epsW) .^ alpha;          % monotone, finite, >0
adaptiveFactor = min(max(confidence, minFactor), maxFactor);   % bounded

illuminated = count(:) > 0;
if any(illuminated)
    mf = mean(adaptiveFactor(illuminated));
    if mf > eps
        adaptiveFactor = adaptiveFactor / mf;      % unit mean -> no global scale distortion
    end
end

weight = geometryWeight(:) .* adaptiveFactor;
weight(~illuminated) = 0;                           % no measurement -> no weight
weight(~isfinite(weight)) = 0;                      % belt-and-braces

d.signalConfidence = signalNorm;
d.adaptiveFactor   = adaptiveFactor;
end
