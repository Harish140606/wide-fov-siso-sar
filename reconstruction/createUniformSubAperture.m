function subAp = createUniformSubAperture(aperture)
% CREATEUNIFORMSUBAPERTURE  Eq. (4)-(7): builds a uniformly-sampled
% sub-aperture from the raw, non-uniformly-dense curved-array
% measurements, by nearest-neighbor selection of actually-measured
% samples closest to an "ideal" uniform grid.
%
%   aperture : struct produced by generateCurvedArray.m (must contain
%              phi_i, theta_j, N, M, R0, pos, beamDir, and the
%              column-major (i,j) linear-index layout used there)
%
%   subAp.idx      : indices into the RAW aperture arrays selected for
%                     the uniform sub-aperture
%   subAp.theta    : selected theta (rad)
%   subAp.phi      : selected phi (rad)
%   subAp.pos      : selected Cartesian positions
%   subAp.beamDir  : selected boresight unit vectors
%   subAp.Mideal   : N x 1, ideal azimuth sample count per elevation (Eq. 5)
%   subAp.d        : reference spatial sampling distance (Eq. 4)

N  = aperture.N;
M  = aperture.M;
R0 = aperture.R0;
phi_i   = aperture.phi_i;
theta_j = aperture.theta_j;

phi_N = phi_i(end);   % phi at i = N, i.e. the array's central angle
if abs(sin(phi_N)) < 1e-12
    error('createUniformSubAperture:zeroPhiN', ...
        'phi_N (max elevation / central angle) must be nonzero (Eq. 4 requires phi_N ~= 0).');
end

d = 2*pi*R0*sin(phi_N) / N;                       % Eq. (4)

Mideal = round(N * sin(phi_i) / sin(phi_N));      % Eq. (5)
Mideal = max(Mideal, 1);                          % phi_i = 0 collapses to a single physical
                                                   % point (the pole) -- keep exactly one sample
                                                   % there instead of dividing by zero azimuth bins

selMask = false(N, M);   % matches the column-major (i,j) layout used by generateCurvedArray

for i = 1:N
    Mi = Mideal(i);
    thetaIdeal = (0:Mi-1)' * (2*pi/Mi);           % Eq. (6)

    % Circular distance from every ideal angle to every actual sampled
    % azimuth angle, then nearest-neighbor selection (Eq. 7).
    diffRaw = thetaIdeal - theta_j';              % Mi x M
    circDist = abs(mod(diffRaw + pi, 2*pi) - pi);
    [~, nearestJ] = min(circDist, [], 2);         % Mi x 1
    nearestJ = unique(nearestJ);

    selMask(i, nearestJ) = true;
end

idxRaw = find(selMask(:));   % column-major linear index, consistent with generateCurvedArray's (:) layout

subAp.idx     = idxRaw;
subAp.theta   = aperture.theta(idxRaw);
subAp.phi     = aperture.phi(idxRaw);
subAp.pos     = aperture.pos(idxRaw, :);
subAp.beamDir = aperture.beamDir(idxRaw, :);
subAp.Mideal  = Mideal;
subAp.d       = d;
subAp.N       = N;

end
