function s = generatePointTargetEcho(aperturePos, targetPos, amplitude, k, beamDir, HPBW_deg)
% GENERATEPOINTTARGETECHO  Single-scatterer special case of Eq. (1),
% kept as its own function because the IPSF / resolution experiments
% (Fig. 5, Section IV) are defined for exactly one point scatterer and
% this makes that correspondence explicit and easy to unit-test.
%
%   aperturePos : Ns x 3
%   targetPos   : 1 x 3
%   amplitude   : scalar complex scattering coefficient (default 1)
%   k           : Nf x 1
%   beamDir, HPBW_deg : optional antenna directivity taper, see
%                 generateSceneEcho.m / antennaGain.m
%   s           : Ns x Nf

if nargin < 3 || isempty(amplitude)
    amplitude = 1;
end
if nargin < 5
    beamDir = [];
    HPBW_deg = [];
end

targets.pos = targetPos;
targets.amp = amplitude;
s = generateSceneEcho(aperturePos, targets, k, beamDir, HPBW_deg);
end
