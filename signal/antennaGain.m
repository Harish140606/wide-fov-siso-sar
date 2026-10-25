function gain = antennaGain(aperturePos, beamDir, targetPos, HPBW_deg)
% ANTENNAGAIN  Physical single-antenna directivity taper, used when
% generating echo data (see generateSceneEcho.m). NOT one of the
% paper's numbered equations -- it is our explicit implementation of
% the physical antenna directivity implied by the paper stating both
% systems use antennas with a 60 deg HPBW (Section IV), which the
% paper's own Eq. (1)/(2) echo model does not otherwise encode.
%
% A cosine-power amplitude pattern g(theta) = max(cos(theta),0)^p is
% used, with theta the angle between an element's boresight and the
% direction to the target, and p chosen so g(HPBW/2) corresponds to a
% -3 dB amplitude point:
%   p = log(10^(-3/20)) / log(cos(HPBW_deg/2 in rad))
%
%   aperturePos, beamDir : Ns x 3
%   targetPos             : 1 x 3
%   HPBW_deg              : scalar
%   gain                  : Ns x 1, in [0,1]

dirVec = targetPos - aperturePos;              % Ns x 3
nrm = sqrt(sum(dirVec.^2, 2));
nrm(nrm < eps) = eps;
cosAng = sum(dirVec .* beamDir, 2) ./ nrm;
cosAng = min(1, max(-1, cosAng));

halfBW = deg2rad(HPBW_deg) / 2;
p = log(10^(-3/20)) / log(cos(halfBW));

gain = max(cosAng, 0) .^ p;
end
