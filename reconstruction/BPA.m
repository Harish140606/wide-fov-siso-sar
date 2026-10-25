function img = BPA(echo, aperturePos, k, voxelPos, cfg, beamDir)
% BPA  Classical back-projection algorithm, Eq. (8), evaluated over the
% RAW (non-uniformized, non-gated) aperture. This is the baseline
% algorithm the paper shows to be inadequate for the curved array
% because of its non-uniform spatial sampling density (Section I,
% Section III intro; Fig. 3).
%
%   echo        : Ns x Nf complex echo data
%   aperturePos : Ns x 3
%   k           : Nf x 1
%   voxelPos    : Nv x 3 (flattened list of image-grid points)
%   cfg         : configuration struct
%   beamDir     : Ns x 3, optional (only needed if you also want to
%                 apply HPBW gating with classical BPA for comparison;
%                 not used in the default un-gated classical BPA path)
%
%   img : Nv x 1 complex reconstructed scattering field

if nargin < 6 || isempty(beamDir)
    beamDir = zeros(size(aperturePos));  % unused when gateOn = false
end

img = backProjectionCore(aperturePos, beamDir, echo, k, voxelPos, cfg, false);
end
