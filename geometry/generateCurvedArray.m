function aperture = generateCurvedArray(cfg)
% GENERATECURVEDARRAY  Eq. (3): spatial sampling coordinates of the
% rotated curved linear array (the "raw", non-uniformized aperture).
%
% The curved SISO array has N fixed elements distributed over the arc
% (elevation angles phi_i, i = 1..N), spanning central angle phi_c.
% The whole array mechanically rotates about the z-axis through M
% azimuth steps covering 360 deg. Every element measures at every
% rotation step, so the raw aperture has N*M spatial samples.
%
%   phi_i   = (i-1)*phi_c/(N-1),  i = 1..N     (elevation, Eq. 3)
%   theta_j = (j-1)*360/M,        j = 1..M     (azimuth,   Eq. 3)
%
% Each antenna's beam points perpendicular to the local tangent of the
% arc, i.e. radially outward from the origin -- since the element sits
% at radius R0 along direction (theta,phi), its own position vector
% already IS that outward unit direction (scaled by R0). This is used
% later as r1 in Eq. (10) and as the HPBW boresight direction.
%
% aperture fields:
%   theta, phi   : M*N x 1 (radians), sampling angles per point
%   x,y,z        : M*N x 1, Cartesian sampling positions
%   pos          : M*N x 3
%   beamDir      : M*N x 3, unit boresight vectors (== pos/R0)
%   elevIdx      : M*N x 1, index i (1..N) of each sample
%   aziIdx       : M*N x 1, index j (1..M) of each sample
%   phi_i, theta_j : unique elevation / azimuth angle lists (rad)
%   R0, N, M

R0    = cfg.curved.R0;
phi_c = deg2rad(cfg.curved.phi_c_deg);
N     = cfg.curved.N;
M     = cfg.curved.M;

i = (1:N)';
phi_i = (i - 1) .* phi_c ./ (N - 1);       % Eq. (3)

j = (1:M)';
theta_j = (j - 1) .* (2*pi) ./ M;          % Eq. (3), 360 deg = 2*pi

[II, JJ] = ndgrid(1:N, 1:M);               % every element x every rotation step
elevIdx = II(:);
aziIdx  = JJ(:);
Phi   = phi_i(elevIdx);
Theta = theta_j(aziIdx);

pos = sphericalToCartesian(R0, Theta, Phi);

aperture.theta   = Theta;
aperture.phi     = Phi;
aperture.x       = pos(:,1);
aperture.y       = pos(:,2);
aperture.z       = pos(:,3);
aperture.pos     = pos;
aperture.beamDir = pos ./ R0;              % radial unit vector, Eq. (10)-consistent
aperture.elevIdx = elevIdx;
aperture.aziIdx  = aziIdx;
aperture.phi_i   = phi_i;
aperture.theta_j = theta_j;
aperture.R0 = R0;
aperture.N  = N;
aperture.M  = M;

end
