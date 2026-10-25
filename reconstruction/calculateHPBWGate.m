function delta = calculateHPBWGate(aperturePos, beamDir, voxelPos, HPBW_deg)
% CALCULATEHPBWGATE  Eq. (10) and Eq. (12): beam-angle gating mask.
%
% HPBW_deg may be:
%   * a scalar                 -> fixed nominal gate (original behaviour), or
%   * an Ns x 1 vector         -> OUR EXTENSION: per-sample effective HPBW
%                                 from computeAdaptiveHPBW.m. Passing a
%                                 constant vector reproduces the scalar case.
%
%   r1 = aperturePos(s,:)          (Eq. 10; for the curved array this
%                                    already equals R0*(sin(phi)cos(theta),
%                                    sin(phi)sin(theta), cos(phi)) since
%                                    the element sits at radius R0 along
%                                    its own boresight direction)
%   r2 = voxelPos - r1             (Eq. 10)
%   alpha(r1,r2) = angle(beamDir, r2)
%   delta = 1 if alpha < HPBW/2, else 0                      (Eq. 12)
%
%   aperturePos, beamDir : Ns x 3 (a chunk of aperture samples and
%                           their boresight unit vectors)
%   voxelPos             : Nv x 3
%   HPBW_deg             : scalar, half-power beamwidth in degrees
%
%   delta : Ns x Nv logical gating matrix

Ns = size(aperturePos, 1);
Nv = size(voxelPos, 1);
halfBW = deg2rad(HPBW_deg(:)) / 2;            % scalar or Ns x 1
if isscalar(halfBW)
    halfBW = repmat(halfBW, Ns, 1);
end

delta = false(Ns, Nv);
for s = 1:Ns
    r2  = voxelPos - aperturePos(s, :);       % Eq. (10), Nv x 3
    nr2 = sqrt(sum(r2.^2, 2));
    nr2(nr2 < eps) = eps;
    cosA = (r2 * beamDir(s, :)') ./ nr2;
    cosA = min(1, max(-1, cosA));
    ang  = acos(cosA);
    delta(s, :) = ang < halfBW(s);            % Eq. (12)
end

end
