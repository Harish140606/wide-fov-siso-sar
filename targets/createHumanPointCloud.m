function human = createHumanPointCloud(yOffset, zDepth, amplitude)
% CREATEHUMANPOINTCLOUD  Simplified, programmatically generated
% "walking human" point-scatterer model (stands in for Fig. 7 of the
% paper). This is NOT a scan of a real body or an external dataset --
% it is a small set of spheres/cylinders approximating head, torso,
% pelvis, arms and legs in a walking pose, built purely from geometry.
%
% AXIS CONVENTION (our own choice, documented here and in the README,
% since the paper does not spell out its own body-relative axes):
%   x : body HEIGHT axis (feet ~ -0.85, head ~ +0.85), i.e. vertical
%   y : lateral standing position (this is the paper's "y = 1.1 m,
%       0 m, -0.8 m" placement) plus small shoulder-width variation
%   z : depth / boresight axis (paper's "z = 1 m" standoff) plus small
%       chest-to-back body-thickness variation
%
%   yOffset : lateral placement of the person, m (paper: 1.1, 0, -0.8)
%   zDepth  : standoff depth of the person, m (paper: 1)
%   amplitude : per-scatterer amplitude (default 1)
%
%   human.pos : Nt x 3
%   human.amp : Nt x 1

if nargin < 3 || isempty(amplitude)
    amplitude = 1;
end

pts = [];

% Point counts kept modest (~600 total) so the human-cloud experiment's
% coherent echo sum (looped per scatterer, see generateSceneEcho.m)
% runs in a reasonable time; increase if higher visual fidelity is
% needed and runtime is not a concern (see USER_GUIDE.md).

% Head
pts = [pts; genSphere([0.78 0 0], 0.09, 30)];
% Torso (chest to hips)
pts = [pts; genCylinder([0.55 0 0], [0.15 0 0], 0.16, 20, 6)];
% Pelvis
pts = [pts; genSphere([0.05 0 0], 0.13, 15)];
% Left leg forward, right leg back (mid-stride walking pose)
pts = [pts; genCylinder([0.02 0.06 0.01], [-0.42 0.16 0.05], 0.065, 12, 5)];   % L thigh
pts = [pts; genCylinder([-0.42 0.16 0.05], [-0.85 0.20 -0.03], 0.05, 12, 5)];  % L shin/foot
pts = [pts; genCylinder([0.02 -0.06 -0.01], [-0.40 -0.14 -0.10], 0.065, 12, 5)]; % R thigh
pts = [pts; genCylinder([-0.40 -0.14 -0.10], [-0.85 -0.10 -0.18], 0.05, 12, 5)]; % R shin/foot
% Left arm back, right arm forward (opposite swing to legs)
pts = [pts; genCylinder([0.50 0.19 0.00], [0.20 0.24 -0.14], 0.045, 10, 5)];   % L upper arm
pts = [pts; genCylinder([0.20 0.24 -0.14], [-0.05 0.26 -0.22], 0.04, 10, 5)];  % L forearm
pts = [pts; genCylinder([0.50 -0.19 0.00], [0.22 -0.22 0.16], 0.045, 10, 5)];  % R upper arm
pts = [pts; genCylinder([0.22 -0.22 0.16], [-0.03 -0.24 0.24], 0.04, 10, 5)];  % R forearm

% Apply lateral / depth placement (height axis x is left untouched)
pts(:,2) = pts(:,2) + yOffset;
pts(:,3) = pts(:,3) + zDepth;

human.pos = pts;
human.amp = amplitude * ones(size(pts,1), 1);

end

% ----------------------------------------------------------------
function pts = genSphere(center, radius, nPts)
% Fibonacci-sphere sampling for an even point distribution on a shell.
idx = (0:nPts-1)' + 0.5;
phi = acos(1 - 2*idx/nPts);
goldenAngle = pi*(1 + sqrt(5));
theta = mod(idx*goldenAngle, 2*pi);
x = radius*sin(phi).*cos(theta);
y = radius*sin(phi).*sin(theta);
z = radius*cos(phi);
pts = center + [x, y, z];
end

% ----------------------------------------------------------------
function pts = genCylinder(p1, p2, radius, nAxial, nRadial)
% Point-shell cylinder from p1 to p2 with the given cross-section radius,
% used as a simple limb/torso primitive.
p1 = p1(:)'; p2 = p2(:)';
axisVec = p2 - p1;
L = norm(axisVec);
if L < eps
    pts = p1;
    return;
end
axisHat = axisVec / L;

% Orthonormal basis perpendicular to the segment axis
tmp = [1 0 0];
if abs(dot(tmp, axisHat)) > 0.9
    tmp = [0 1 0];
end
e1 = cross(axisHat, tmp); e1 = e1 / norm(e1);
e2 = cross(axisHat, e1);

t = linspace(0, 1, nAxial)';
ang = linspace(0, 2*pi, nRadial+1)'; ang(end) = [];

pts = zeros(nAxial*numel(ang), 3);
n = 0;
for a = 1:numel(ang)
    ring = p1 + t*axisVec + radius*cos(ang(a))*e1 + radius*sin(ang(a))*e2;
    pts(n+1:n+nAxial, :) = ring;
    n = n + nAxial;
end
end
