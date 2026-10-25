function [effHPBW_deg, summary] = computeAdaptiveHPBW(subAp, cfg)
% COMPUTEADAPTIVEHPBW  OUR PROPOSED EXTENSION (NOT from the paper).
%
% Produces a PER-SAMPLE effective half-power beamwidth (deg) for the
% Eq.(12) gate, from sampling geometry only -- no signal feedback, so
% there is no circular dependency on the reconstruction and a single
% pass is sufficient.
%
% Rationale: on the curved array the central elevation rings (small phi)
% are heavily over-sampled and produce the central circular clutter the
% paper describes; the peripheral rings (phi near phi_c) are weakly
% sampled and lose coverage. We therefore gate the central rings
% slightly TIGHTER (reject more off-axis clutter) and the peripheral
% rings slightly WIDER (retain more of their sparse contributions),
% while keeping every value a physically meaningful beamwidth bounded by
% the configured limits.
%
%   normElev(s)     = phi(s) / max(phi)                    in [0, 1]
%   effHPBW_deg(s)  = HPBW_deg + adaptStrength*HPBW_deg*(normElev(s) - 0.5)
%   effHPBW_deg(s)  = clip(effHPBW_deg(s), minHPBW_deg, maxHPBW_deg)
%
% adaptStrength = 0  ->  effHPBW_deg == cfg.HPBW_deg for every sample,
% i.e. identical to the fixed nominal gate.
%
%   subAp : uniform sub-aperture struct (needs .phi, Ns x 1, radians)
%   cfg   : configuration struct (cfg.HPBW_deg, cfg.gating.*)
%
%   effHPBW_deg : Ns x 1, effective per-sample HPBW in degrees
%   summary     : struct (mode, nominal_deg, min_deg, max_deg, mean_deg)

phi = subAp.phi(:);
Ns  = numel(phi);

nominal = cfg.HPBW_deg;
strength = cfg.gating.adaptStrength;
loLim = cfg.gating.minHPBW_deg;
hiLim = cfg.gating.maxHPBW_deg;

phiMax = max(phi);
if phiMax < eps
    normElev = zeros(Ns, 1);
else
    normElev = phi / phiMax;
end

effHPBW_deg = nominal + strength * nominal * (normElev - 0.5);
effHPBW_deg = min(max(effHPBW_deg, loLim), hiLim);   % hard physical bounds

summary = struct('mode', 'adaptive', ...
                 'nominal_deg', nominal, ...
                 'min_deg', min(effHPBW_deg), ...
                 'max_deg', max(effHPBW_deg), ...
                 'mean_deg', mean(effHPBW_deg));
end
