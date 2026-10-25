function frac = computeArtifactEnergy(imgMag, trueMask)
% COMPUTEARTIFACTENERGY  ADDITIONAL metric (NOT from the paper).
% Fraction of total image energy that falls OUTSIDE the known true
% target support (e.g. outside the point-array grid locations) --
% a proxy for the "circular artifact / sidelobe" energy the paper
% describes qualitatively. Lower = less artifact energy.
%
%   imgMag   : image magnitude (any shape)
%   trueMask : logical, same shape as imgMag, true at genuine target
%              support locations (a small neighborhood around each
%              known scatterer)

e = imgMag(:).^2;
total = sum(e);
if total < eps
    frac = NaN;
    return;
end
artifact = sum(e(~trueMask(:)));
frac = artifact / total;
end
