function s = generateSceneEcho(aperturePos, targets, k, beamDir, HPBW_deg)
% GENERATESCENEECHO  Discrete-scatterer implementation of Eq. (1).
%
% The continuous integral
%   s(theta,phi,k) = INT sigma(x',y',z') * exp(-2j*k*RR) dx' dy' dz'
% is realized as a coherent sum over a finite set of point scatterers
% (the standard discretization used for numerical scattering
% simulation under the first-order Born approximation, Eq. (1)-(2)).
%
%   aperturePos : Ns x 3, spatial sampling positions (R0,theta,phi) in
%                 Cartesian form
%   targets.pos : Nt x 3, scatterer positions
%   targets.amp : Nt x 1, scatterer complex amplitudes (scattering
%                 coefficients sigma)
%   k           : Nf x 1, wavenumber grid
%   beamDir, HPBW_deg : OPTIONAL. If given, a per-element antenna
%                 directivity taper is applied to the echo amplitude,
%                 attenuating a scatterer's contribution as its
%                 direction moves away from that element's boresight.
%                 This models the physical fact -- stated in Section IV
%                 for BOTH the curved and the planar comparison system
%                 ("...60 deg HPBW") -- that a real antenna element has
%                 finite directivity; without it, an ideal isotropic
%                 element would "see" every target from every angle,
%                 which is why the planar system's degraded/failed
%                 imaging at large angles (Fig. 5) does not otherwise
%                 emerge from the geometry alone. A standard cosine-power
%                 pattern is used, calibrated so the amplitude is -3 dB
%                 at +/- HPBW/2 off boresight (see antennaGain.m). This
%                 is separate from, and applied in addition to, the
%                 Eq. (12) HPBW gate used only by the proposed
%                 reconstruction algorithm.
%
%   s : Ns x Nf complex echo matrix, s(sample, freq)

if nargin < 4
    beamDir = [];
    HPBW_deg = [];
end

Ns = size(aperturePos, 1);
Nf = numel(k);
Nt = size(targets.pos, 1);

s = zeros(Ns, Nf);
kRow = k(:)';   % 1 x Nf

for t = 1:Nt
    RR = calculateRange(aperturePos, targets.pos(t,:));  % Ns x 1, Eq. (2)
    if ~isempty(beamDir)
        gain = antennaGain(aperturePos, beamDir, targets.pos(t,:), HPBW_deg);  % Ns x 1
    else
        gain = 1;
    end
    s = s + targets.amp(t) .* gain .* exp(-2i * RR * kRow);       % Eq. (1)
end
end
