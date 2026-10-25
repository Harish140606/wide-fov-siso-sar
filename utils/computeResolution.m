function width = computeResolution(axisVals, profileMag)
% COMPUTERESOLUTION  Metric directly corresponding to the paper's
% methodology (Section IV, Fig. 5 caption): the -3 dB width of an IPSF
% line slice along one axis.
%
%   axisVals   : 1 x P or P x 1, coordinate of each sample along the cut
%   profileMag : magnitude of the reconstructed field along that cut
%   width      : -3 dB width, same units as axisVals. NaN if the
%                profile never drops 3 dB from its peak within the
%                provided window (peak too close to the window edge).

axisVals = axisVals(:);
profileMag = profileMag(:);

pk = max(profileMag);
if pk < eps
    width = NaN;
    return;
end
profileDB = 20*log10(profileMag/pk + eps);

[~, peakIdx] = max(profileMag);

leftIdx = find(profileDB(1:peakIdx) <= -3, 1, 'last');
rightRel = find(profileDB(peakIdx:end) <= -3, 1, 'first');

if isempty(leftIdx) || isempty(rightRel)
    width = NaN;
    return;
end
rightIdx = rightRel + peakIdx - 1;

xL = interp1(profileDB(leftIdx:leftIdx+1), axisVals(leftIdx:leftIdx+1), -3, 'linear');
xR = interp1(profileDB(rightIdx-1:rightIdx), axisVals(rightIdx-1:rightIdx), -3, 'linear');

width = xR - xL;
end
