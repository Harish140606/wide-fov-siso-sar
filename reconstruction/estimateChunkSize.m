function [chunkSize, estMB] = estimateChunkSize(Nv, Ns, Nf, cfg)
% ESTIMATECHUNKSIZE  OUR EXTENSION. Memory-aware automatic aperture
% chunk-size selection for backProjectionCore.m.
%
% The per-chunk working set of the FFT path is dominated by:
%   G          : chunkSize x Nfft   complex  (IFFT of the chunk's echo)
%   per-voxel  : a handful of  Nv x 1  temporaries (RR, tau, delta, ...)
%
% We size the chunk so that
%   16 bytes * (chunkSize*Nfft + safety*Nv)  <=  targetMemoryMB
% and clamp the result to [1, Ns] as well as to the configured fixed
% chunkSize is chosen purely from the memory budget. If a system memory
% figure is available (Windows
% `memory`), the target is additionally capped at a fraction of the free
% physical memory; if unavailable we fall back silently to the configured
% cfg.performance.targetMemoryMB, which always exists.
%
%   Nv, Ns, Nf : problem sizes
%   cfg        : configuration struct (uses cfg.recon.zeroPadFactor,
%                cfg.recon.chunkSize, cfg.performance.targetMemoryMB)
%
%   chunkSize : positive integer, 1 <= chunkSize <= Ns
%   estMB     : estimated per-chunk working-set size at that chunkSize

bytesPerComplex = 16;
safetyVoxTemps  = 8;      % ~8 full Nv-length temporaries live at once
overheadFactor  = 1.5;    % copies / MATLAB internal slack

zpf  = cfg.recon.zeroPadFactor;
Nfft = max(zpf * Nf, Nf);

targetMB = cfg.performance.targetMemoryMB;
try
    m = memory;  %#ok<MMEM>  Windows only
    freeMB = m.MemAvailableAllArrays / 1e6;
    targetMB = min(targetMB, 0.25 * freeMB);
catch
    % no `memory` (non-Windows) -> keep configured budget
end
targetBytes = targetMB * 1e6 / overheadFactor;

fixedTerm = bytesPerComplex * safetyVoxTemps * Nv;
avail = targetBytes - fixedTerm;
if avail < bytesPerComplex * Nfft
    chunkSize = 1;                       % extreme case: barely fits one sample
else
    chunkSize = floor(avail / (bytesPerComplex * Nfft));
end

chunkSize = max(1, min(chunkSize, Ns));
estMB = bytesPerComplex * (chunkSize * Nfft + safetyVoxTemps * Nv) * overheadFactor / 1e6;

end
