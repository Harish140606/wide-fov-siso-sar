function pos = sphericalToCartesian(r, theta, phi)
% SPHERICALTOCARTESIAN  Eq. (3) coordinate transform.
%
% pos = sphericalToCartesian(r, theta, phi)
%   r, theta, phi : scalars or equal-size arrays. theta (azimuth) and
%                   phi (elevation) MUST be in radians. r may be a
%                   scalar (curved array radius R0) or an array.
%   pos           : N x 3 matrix of [x y z] Cartesian coordinates,
%                   N = numel(theta).
%
% x = r*sin(phi)*cos(theta)
% y = r*sin(phi)*sin(theta)
% z = r*cos(phi)

theta = theta(:);
phi   = phi(:);
if isscalar(r)
    r = r * ones(size(theta));
else
    r = r(:);
end

x = r .* sin(phi) .* cos(theta);
y = r .* sin(phi) .* sin(theta);
z = r .* cos(phi);

pos = [x, y, z];
end
