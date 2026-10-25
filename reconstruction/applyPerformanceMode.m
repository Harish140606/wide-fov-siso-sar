function [cfg, info] = applyPerformanceMode(cfg, Nv, Ns, Nf)
% APPLYPERFORMANCEMODE  OUR EXTENSION (implementation contribution).
% Resolves the high-level cfg.performance.* section onto the existing
% low-level engine knobs cfg.recon.method / cfg.useGPU / cfg.useParallel
% / cfg.recon.chunkSize that backProjectionCore.m actually consumes.
%
% This function does NOT change any numerical result -- it only chooses
% an evaluation strategy. With the default cfg.performance settings
% (mode='fft_serial', chunkMode='fixed', chunkSize=0) it is a no-op:
% method stays 'fft', GPU/parallel stay off, chunkSize is untouched, so
% the original paper behaviour is reproduced exactly.
%
% If cfg has no .performance field at all (e.g. an externally hand-built
% cfg), this function leaves cfg completely untouched -- old callers are
% unaffected.
%
%   cfg            : configuration struct
%   Nv, Ns, Nf     : #voxels, #aperture samples, #frequencies (only used
%                    for chunkMode='auto'; may be omitted)
%
%   cfg  : same struct with cfg.recon.method / cfg.useGPU /
%          cfg.useParallel / cfg.recon.chunkSize resolved
%   info : struct describing the decision (requestedMode, resolvedMode,
%          chunkMode, chunkSize, autoReason)

info = struct('requestedMode', '(none)', 'resolvedMode', cfg.recon.method, ...
              'chunkMode', 'fixed', 'chunkSize', cfg.recon.chunkSize, 'autoReason', '');

if ~isfield(cfg, 'performance') || ~isstruct(cfg.performance)
    return;   % nothing to do -- preserve legacy behaviour
end
p = cfg.performance;

% ---- 1) hardware / method mode --------------------------------------
requested = lower(strtrim(p.mode));
info.requestedMode = requested;

resolved = requested;
if strcmp(resolved, 'fft_auto')
    [resolved, info.autoReason] = pickAutoMode();
end

switch resolved
    case 'direct'
        cfg.recon.method = 'direct'; cfg.useGPU = false; cfg.useParallel = false;
    case 'fft_serial'
        cfg.recon.method = 'fft';    cfg.useGPU = false; cfg.useParallel = false;
    case 'fft_parallel'
        cfg.recon.method = 'fft';    cfg.useGPU = false; cfg.useParallel = true;
    case 'fft_gpu'
        cfg.recon.method = 'fft';    cfg.useGPU = true;  cfg.useParallel = false;
    otherwise
        error('applyPerformanceMode:badMode', ...
            'cfg.performance.mode must be one of: direct | fft_serial | fft_parallel | fft_gpu | fft_auto (got ''%s'').', p.mode);
end
info.resolvedMode = resolved;

% ---- 2) chunk sizing ----------------------------------------------
info.chunkMode = lower(strtrim(p.chunkMode));
if strcmp(info.chunkMode, 'auto') && nargin >= 4
    cfg.recon.chunkSize = estimateChunkSize(Nv, Ns, Nf, cfg);
elseif isfield(p, 'chunkSize') && ~isempty(p.chunkSize) && p.chunkSize > 0
    cfg.recon.chunkSize = round(p.chunkSize);
end
cfg.recon.chunkSize = max(1, round(cfg.recon.chunkSize));
info.chunkSize = cfg.recon.chunkSize;

end

% ------------------------------------------------------------------
function [mode, reason] = pickAutoMode()
% Choose the fastest SAFE available FFT mode from the hardware/toolboxes
% actually present. Order of preference: GPU > parallel CPU > serial.
mode = 'fft_serial';
reason = 'no GPU and no Parallel Computing Toolbox detected';

hasPCT = ~isempty(ver('parallel'));

gpuOK = false;
if hasPCT
    try
        gpuOK = gpuDeviceCount() > 0;
    catch
        gpuOK = false;
    end
end

if gpuOK
    mode = 'fft_gpu';
    reason = 'supported GPU detected';
elseif hasPCT
    mode = 'fft_parallel';
    reason = 'Parallel Computing Toolbox detected, no GPU';
end
end
