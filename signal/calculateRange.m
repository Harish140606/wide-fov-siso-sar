function RR = calculateRange(P1, P2)
% CALCULATERANGE  Eq. (2): Euclidean distance between two point sets.
%
%   RR = calculateRange(P1, P2)
%   P1 : N1 x 3, P2 : N2 x 3  (Cartesian [x y z])
%   RR : N1 x N2 matrix of pairwise Euclidean distances,
%        RR(a,b) = || P1(a,:) - P2(b,:) ||
%
% This is the generic pairwise-distance evaluation of
%   RR = sqrt((x-xs)^2 + (y-ys)^2 + (z-zs)^2)
% used both for target-to-sampling-point ranges (echo generation) and
% for voxel-to-sampling-point ranges (reconstruction). Which set plays
% "target" vs "sampling point" is determined by the caller.

dx = P1(:,1) - P2(:,1)';
dy = P1(:,2) - P2(:,2)';
dz = P1(:,3) - P2(:,3)';
RR = sqrt(dx.^2 + dy.^2 + dz.^2);
end
