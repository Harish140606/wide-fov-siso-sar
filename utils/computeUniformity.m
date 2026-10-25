function u = computeUniformity(peakVals)
% COMPUTEUNIFORMITY  ADDITIONAL metric (NOT from the paper). Quantifies
% how evenly a set of point-scatterer reconstruction peaks (e.g. the
% planar point-array experiment) are recovered. Lower coefficientOfVariation
% = more uniform reconstruction across the scene -- used to give a
% number to the paper's qualitative claim of "more uniform" images.

peakVals = abs(peakVals(:));
u.mean = mean(peakVals);
u.std  = std(peakVals);
if u.mean > eps
    u.coefficientOfVariation = u.std / u.mean;
else
    u.coefficientOfVariation = NaN;
end
end
