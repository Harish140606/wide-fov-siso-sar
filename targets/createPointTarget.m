function target = createPointTarget(r, theta_deg, phi_deg, amplitude)
% CREATEPOINTTARGET  Single point scatterer at spherical coordinates
% (r, theta, phi), converted to Cartesian via the same transform used
% for the array geometry (sphericalToCartesian, Eq. 3 form).
%
%   target.pos : 1 x 3
%   target.amp : scalar

if nargin < 4 || isempty(amplitude)
    amplitude = 1;
end

theta = deg2rad(theta_deg);
phi   = deg2rad(phi_deg);

target.pos = sphericalToCartesian(r, theta, phi);
target.amp = amplitude;
end
