function aperture = generatePlanarArray(cfg)
% GENERATEPLANARARRAY  Comparison "classic" planar-scanning SISO-SAR
% (Fig. 1(a) of the paper): a 1-D linear array of length L placed along
% x, mechanically scanned along y, forming a 2-D planar aperture at
% z = 0, boresight fixed along +z (front-looking).
%
% This system is NOT governed by Eq. (3) (that equation is specific to
% the curved/rotated array); it is our direct, literal implementation
% of the "planar scanning" baseline geometry described in Section I and
% used for comparison in Section IV.
%
% Parameters (Section IV): array length 0.6 m, 353 Tx/Rx elements,
% 353 scan positions.

L  = cfg.planar.length;
Ne = cfg.planar.numElements;
Ns = cfg.planar.numScan;

xElem = linspace(-L/2, L/2, Ne)';
yScan = linspace(-L/2, L/2, Ns)';

[XX, YY] = ndgrid(xElem, yScan);
x = XX(:);
y = YY(:);
z = zeros(numel(x), 1);
pos = [x, y, z];

aperture.x = x;
aperture.y = y;
aperture.z = z;
aperture.pos = pos;
aperture.beamDir = repmat([0 0 1], numel(x), 1);  % fixed front-looking boresight
aperture.N = Ne;
aperture.M = Ns;

end
