function weight = calculateSpatialWeights(aperturePos, beamDir, voxelPos, HPBW_deg, chunkSize)
% CALCULATESPATIALWEIGHTS  Eq. (9): spatial-grid weight.
%
%   w(x,y,z) = 1 / |{(theta,phi) : alpha(r1,r2) < Btheta/2}|
%
% i.e. the reciprocal of the number of (uniform sub-aperture) samples
% that illuminate a given voxel within the antenna half-power
% beamwidth. Voxels seen by many samples (typically the over-sampled
% center of the hemispherical aperture) are down-weighted; voxels seen
% by few samples (the periphery) are boosted -- this is the
% "weak-illumination area enhancement" described in Section III/I.
%
% This function depends only on sampling geometry and the antenna
% beamwidth (not on the measured echo or frequency), so -- exactly as
% stated in the paper ("weights are computed only once per imaging
% scene") -- it only needs to be evaluated once per (aperture, voxel
% grid) pair and reused for every frequency/target.
%
%   aperturePos, beamDir : Ns x 3 (the UNIFORM sub-aperture, i.e. the
%                           output of createUniformSubAperture.m)
%   voxelPos              : Nv x 3
%   HPBW_deg               : scalar (fixed nominal gate, original), OR
%                            an Ns x 1 vector of per-sample effective
%                            HPBW values (OUR EXTENSION, from
%                            computeAdaptiveHPBW.m). The vector case must
%                            use the SAME per-sample HPBW that the gated
%                            back-projection uses, so that `count` here
%                            matches the illumination count there.
%   chunkSize               : aperture samples processed per chunk
%
%   weight : Nv x 1. Voxels illuminated by zero samples get weight 0
%            (not Inf/NaN) -- there is no measurement to weight there.

Ns = size(aperturePos, 1);
Nv = size(voxelPos, 1);

perSample = ~isscalar(HPBW_deg);

count = zeros(Nv, 1);
for c0 = 1:chunkSize:Ns
    idx = c0:min(c0 + chunkSize - 1, Ns);
    if perSample
        hb = HPBW_deg(idx);
    else
        hb = HPBW_deg;
    end
    delta = calculateHPBWGate(aperturePos(idx, :), beamDir(idx, :), voxelPos, hb); % Eq. (12)
    count = count + sum(delta, 1)';
end

weight = zeros(Nv, 1);
illuminated = count > 0;
weight(illuminated) = 1 ./ count(illuminated);   % Eq. (9)

end
