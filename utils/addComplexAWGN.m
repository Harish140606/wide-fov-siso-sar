function [noisy, noisePower] = addComplexAWGN(echo, SNR_dB, seed)
% ADDCOMPLEXAWGN  OUR EXTENSION. Adds circularly-symmetric complex
% Gaussian noise to an Ns x Nf echo matrix at a specified SNR.
%
% SNR is defined relative to the MEAN power of the (finite, nonzero)
% echo samples:
%   signalPower = mean(|echo|.^2)   over all entries
%   noisePower  = signalPower / 10^(SNR_dB/10)
%   noise       = sqrt(noisePower/2) * (randn + 1j*randn)
%
%   echo    : Ns x Nf complex
%   SNR_dB  : scalar; Inf (or []) -> returns echo unchanged (no noise)
%   seed    : optional RNG seed for reproducibility ([] -> do not reseed)
%
%   noisy      : Ns x Nf complex
%   noisePower : the per-sample noise variance actually used (0 if none)

if nargin < 2 || isempty(SNR_dB) || isinf(SNR_dB)
    noisy = echo;
    noisePower = 0;
    return;
end

if nargin >= 3 && ~isempty(seed)
    rng(seed);
end

sig = echo(:);
sig = sig(isfinite(sig));
signalPower = mean(abs(sig).^2);
if ~isfinite(signalPower) || signalPower <= 0
    noisy = echo;
    noisePower = 0;
    return;
end

noisePower = signalPower / (10^(SNR_dB/10));
n = sqrt(noisePower/2) * (randn(size(echo)) + 1j*randn(size(echo)));
noisy = echo + n;
end
