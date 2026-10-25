function targets = createPointArray(z, halfWidth, spacing, amplitude)
% CREATEPOINTARRAY  Planar grid of point scatterers at fixed depth z
% (reproduces the "planar point array at z = 1 m" scene used in
% Fig. 3(b) / Fig. 6 of the paper).
%
%   z         : depth of the plane, m (paper uses 1 m)
%   halfWidth : plane half-extent in x and y, m
%   spacing   : scatterer spacing, m
%   amplitude : scalar amplitude applied to every scatterer (default 1)
%
%   targets.pos : Nt x 3
%   targets.amp : Nt x 1

if nargin < 4 || isempty(amplitude)
    amplitude = 1;
end

xs = -halfWidth:spacing:halfWidth;
ys = -halfWidth:spacing:halfWidth;
[XX, YY] = ndgrid(xs, ys);

pos = [XX(:), YY(:), z*ones(numel(XX), 1)];
targets.pos = pos;
targets.amp = amplitude * ones(size(pos,1), 1);
end
