function dB = toDB(img, dynRangeDB)
% TODB  Numerically-safe conversion to normalized dB magnitude.
%   dB = 20*log10(|img| / (max(|img|) + eps) + eps)
% Optionally floors the output at -dynRangeDB (matches the paper's
% figures, which display a fixed dB floor, e.g. -15 dB).

mag = abs(img);
m = max(mag(:));
if m < eps
    m = eps;
end
dB = 20*log10(mag/m + eps);

if nargin >= 2 && ~isempty(dynRangeDB)
    dB = max(dB, -dynRangeDB);
end
end
