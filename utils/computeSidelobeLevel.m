function pslDB = computeSidelobeLevel(axisVals, profileMag, mainlobeHalfWidth)
% COMPUTESIDELOBELEVEL  ADDITIONAL metric (NOT reported numerically by
% the paper -- the paper only shows sidelobes qualitatively in its
% images). Peak sidelobe level relative to the main lobe peak, in dB,
% excluding a window of +/- mainlobeHalfWidth around the peak.

axisVals = axisVals(:);
profileMag = profileMag(:);

pk = max(profileMag);
if pk < eps
    pslDB = NaN;
    return;
end
[~, peakIdx] = max(profileMag);
peakPos = axisVals(peakIdx);

outside = abs(axisVals - peakPos) > mainlobeHalfWidth;
if ~any(outside)
    pslDB = -Inf;
    return;
end

pslDB = 20*log10(max(profileMag(outside))/pk + eps);
end
