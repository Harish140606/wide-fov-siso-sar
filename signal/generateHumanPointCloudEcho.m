function s = generateHumanPointCloudEcho(aperturePos, humanCloud, k, beamDir, HPBW_deg)
% GENERATEHUMANPOINTCLOUDECHO  Thin wrapper applying Eq. (1) (via
% generateSceneEcho) to a walking-human point-cloud target set
% (Fig. 7 / Fig. 8 experiment). Kept as a separate named function per
% the project's module map even though the underlying model is
% identical to generateSceneEcho.
%
%   aperturePos : Ns x 3
%   humanCloud  : struct with fields .pos (Nt x 3), .amp (Nt x 1)
%   k           : Nf x 1
%   beamDir, HPBW_deg : optional antenna directivity taper, see
%                 generateSceneEcho.m / antennaGain.m
%   s           : Ns x Nf

if nargin < 4
    beamDir = [];
    HPBW_deg = [];
end
s = generateSceneEcho(aperturePos, humanCloud, k, beamDir, HPBW_deg);
end
