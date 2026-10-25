function [f, k] = generateFrequencyGrid(cfg)
% GENERATEFREQUENCYGRID  Stepped-frequency sweep used in Eq. (1)/(2).
%   f : numFreq x 1, Hz, linearly spaced f_start..f_end
%   k = 2*pi*f/c : numFreq x 1, wavenumber (as defined in Section II)

f = linspace(cfg.freq.f_start, cfg.freq.f_end, cfg.freq.numFreq)';
k = 2*pi*f / cfg.c;
end
