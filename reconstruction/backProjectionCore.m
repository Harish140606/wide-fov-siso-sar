function [img, count, info] = backProjectionCore(aperturePos, beamDir, echo, k, voxelPos, cfg, gateOn, gateHPBW_deg)
% BACKPROJECTIONCORE  Shared numerical engine implementing the
% back-projection / time-domain-correlation integral of Eq. (8):
%
%   sigma_hat(x,y,z) = INT INT INT s(theta,phi,k) * exp(+2j*k*RR) ...
%                        dtheta dphi dk
%
% realized discretely as a coherent sum over aperture samples and
% frequency:
%
%   img(voxel) = sum_samples [ sum_freq  s(sample,f) * exp(2j*k(f)*RR) ]
%
% Two mathematically-equivalent evaluation methods are provided:
%
%   cfg.recon.method = 'direct'
%       Literal brute-force evaluation of the inner sum over k for
%       every (voxel, sample) pair. Exact, but O(Nvoxel*Nsample*Nfreq)
%       -- only usable for small validation problems.
%
%   cfg.recon.method = 'fft'  (default, used for all real experiments)
%       The inner sum over k, for a FIXED sample, is exactly the
%       inverse discrete Fourier transform of s(sample,:) evaluated at
%       "delay" tau = 2*RR/c (this follows directly from substituting
%       k = 2*pi*f/c into exp(2j*k*RR) = exp(j*2*pi*f*(2R/c)), i.e. a
%       standard Fourier delay/range-compression pair). We therefore
%       IFFT each sample's spectrum once (zero-padded for a fine range
%       grid) and linearly interpolate the resulting range profile at
%       tau = 2*RR(voxel,sample)/c for every voxel, instead of summing
%       401 (or however many) frequency terms explicitly per voxel.
%       validate_project.m checks this equivalence numerically.
%
% If gateOn = true, the Eq. (12) half-power-beamwidth gate delta is
% applied during accumulation (this is what turns classical BPA into
% the gated sum used by the proposed measurement-fusion algorithm);
% count(voxel) then returns the number of aperture samples that
% illuminated that voxel (needed for Eq. (9) weighting). If
% gateOn = false, delta == 1 everywhere (classical BPA, Eq. 8 as-is).
%
% Performance (Part 13): the aperture is processed in chunks of
% cfg.recon.chunkSize samples to bound memory. If cfg.useParallel is
% true, chunks are distributed across a parfor loop (requires the
% Parallel Computing Toolbox to be installed -- NOT merely licensed;
% this is checked and silently falls back to a serial loop with a
% one-time warning if unavailable). If cfg.useGPU is true, each
% chunk's arithmetic is done on the GPU via gpuArray (requires the
% Parallel Computing Toolbox and a supported GPU; same safe fallback).
% Neither flag is required for correctness -- both default to false
% and every numeric result in this project was produced with them off.
%
%   aperturePos : Ns x 3
%   beamDir     : Ns x 3 (unit boresight vectors; only used if gateOn)
%   echo        : Ns x Nf complex
%   k           : Nf x 1
%   voxelPos    : Nv x 3
%   cfg         : configuration struct (see config.m)
%   gateOn      : logical scalar
%   gateHPBW_deg: OPTIONAL (OUR EXTENSION). Full HPBW in degrees used for
%                 the Eq.(12) gate. Scalar (default cfg.HPBW_deg, i.e.
%                 the original fixed gate) or an Ns x 1 vector of
%                 per-sample effective HPBW values (adaptive gating).
%
%   img   : Nv x 1 complex reconstructed field
%   count : Nv x 1 illumination counts (only meaningful if gateOn)
%   info  : struct describing the execution strategy actually used
%           (requestedMode, resolvedMode, executionMode, numWorkers,
%            chunkSize, chunkMode) -- OUR EXTENSION, backward compatible.

Ns = size(aperturePos, 1);
Nv = size(voxelPos, 1);

c  = cfg.c;
f  = k(:) * c / (2*pi);
Nf = numel(f);
df = mean(diff(f));

% OUR EXTENSION: resolve cfg.performance.* onto cfg.recon.* / useGPU /
% useParallel. No-op (and hence original behaviour) under the default
% cfg.performance settings, or when cfg has no .performance field.
[cfg, perfInfo] = applyPerformanceMode(cfg, Nv, Ns, Nf);

if nargin < 8 || isempty(gateHPBW_deg)
    gateHPBW_deg = cfg.HPBW_deg;
end
halfBWvec = deg2rad(gateHPBW_deg(:)) / 2;     % scalar or Ns x 1
if isscalar(halfBWvec)
    halfBWvec = repmat(halfBWvec, Ns, 1);
end

method  = lower(cfg.recon.method);
chunkSz = cfg.recon.chunkSize;

if strcmp(method, 'fft')
    Nfft = max(cfg.recon.zeroPadFactor * Nf, Nf);
    dtau = 1 / (Nfft * df);
elseif ~strcmp(method, 'direct')
    error('backProjectionCore:badMethod', 'cfg.recon.method must be ''fft'' or ''direct''.');
else
    Nfft = 0; dtau = 0;
end

useGPU      = isfield(cfg, 'useGPU')      && cfg.useGPU      && canUseGPU();
useParallel = isfield(cfg, 'useParallel') && cfg.useParallel && canUseParallel();

chunkStarts = 1:chunkSz:Ns;
numChunks = numel(chunkStarts);

if useParallel
    imgParts   = zeros(Nv, numChunks);
    countParts = zeros(Nv, numChunks);
    parfor ci = 1:numChunks
        c0 = chunkStarts(ci); %#ok<PFBNS>
        idx = c0:min(c0 + chunkSz - 1, Ns);
        % each worker writes ONLY its own column -> no race condition;
        % the partial images are summed after the loop.
        [imgParts(:,ci), countParts(:,ci)] = processChunk( ...
            aperturePos(idx,:), beamDir(idx,:), echo(idx,:), f, voxelPos, ...
            method, Nfft, dtau, halfBWvec(idx), gateOn, c, useGPU); %#ok<PFBNS>
    end
    img   = sum(imgParts, 2);
    count = sum(countParts, 2);
else
    img   = zeros(Nv, 1);
    count = zeros(Nv, 1);
    for ci = 1:numChunks
        c0 = chunkStarts(ci);
        idx = c0:min(c0 + chunkSz - 1, Ns);
        [imgPart, countPart] = processChunk( ...
            aperturePos(idx,:), beamDir(idx,:), echo(idx,:), f, voxelPos, ...
            method, Nfft, dtau, halfBWvec(idx), gateOn, c, useGPU);
        img   = img + imgPart;
        count = count + countPart;
    end
end

% ---- execution diagnostics (OUR EXTENSION) ------------------------
if useGPU
    execMode = 'fft_gpu';
elseif useParallel
    execMode = 'fft_parallel';
elseif strcmp(method, 'direct')
    execMode = 'direct';
else
    execMode = 'fft_serial';
end

nWorkers = 1;
if useParallel
    try
        pool = gcp('nocreate');
        if ~isempty(pool), nWorkers = pool.NumWorkers; end
    catch
    end
end

info = perfInfo;
info.executionMode = execMode;      % what actually ran (after fallbacks)
info.numWorkers    = nWorkers;
info.numChunks     = numChunks;
info.chunkSize     = chunkSz;
info.method        = method;

if isfield(cfg, 'performance') && isfield(cfg.performance, 'verbose') && cfg.performance.verbose
    fprintf(['[backProjectionCore] requested=%s  resolved=%s  actual=%s  ' ...
             'workers=%d  chunkMode=%s  chunkSize=%d  chunks=%d\n'], ...
        info.requestedMode, info.resolvedMode, info.executionMode, ...
        info.numWorkers, info.chunkMode, info.chunkSize, info.numChunks);
end

end

% ----------------------------------------------------------------
function [imgPart, countPart] = processChunk(apChunk, beamChunk, echoChunk, f, voxelPos, method, Nfft, dtau, halfBW, gateOn, c, useGPU)
% halfBW : per-sample half-beamwidth (rad) for this chunk, length = #samples in chunk

Nv = size(voxelPos, 1);

if useGPU
    apChunk = gpuArray(apChunk);
    beamChunk = gpuArray(beamChunk);
    echoChunk = gpuArray(echoChunk);
    voxelPos = gpuArray(voxelPos);
    f = gpuArray(f);
    imgPart = gpuArray.zeros(Nv, 1, 'like', 1i);
    countPart = gpuArray.zeros(Nv, 1);
else
    imgPart = zeros(Nv, 1);
    countPart = zeros(Nv, 1);
end

if strcmp(method, 'fft')
    G = ifft(echoChunk, Nfft, 2);
end

for s = 1:size(apChunk, 1)
    r2 = voxelPos - apChunk(s, :);              % Eq. (10), Nv x 3
    RR = sqrt(sum(r2.^2, 2));                    % Eq. (2), Nv x 1

    if gateOn
        nr2  = max(RR, eps);
        cosA = (r2 * beamChunk(s, :)') ./ nr2;
        cosA = min(1, max(-1, cosA));
        delta = acos(cosA) < halfBW(s);           % Eq. (12); halfBW(s) == nominal in fixed mode
    else
        delta = true(Nv, 1);
    end

    switch method
        case 'fft'
            tau  = 2 * RR / c;
            bin  = tau / dtau;
            b0   = floor(bin);
            frac = bin - b0;
            b0m  = mod(b0, Nfft) + 1;
            b1m  = mod(b0 + 1, Nfft) + 1;
            gval = G(s, b0m).' .* (1 - frac) + G(s, b1m).' .* frac;
            carrier = exp(1i * 2 * pi * f(1) * tau);
            contrib = gval .* carrier;
        case 'direct'
            contrib = exp(2i * RR * f(:)' * 2 * pi / c) * echoChunk(s, :).';
    end

    imgPart   = imgPart + delta .* contrib;
    countPart = countPart + delta;
end

if useGPU
    imgPart = gather(imgPart);
    countPart = gather(countPart);
end
end

% ----------------------------------------------------------------
function tf = canUseParallel()
persistent cached
if isempty(cached)
    try
        gcp('nocreate');
        cached = true;
    catch
        warning('backProjectionCore:noParallel', ...
            'cfg.useParallel=true but the Parallel Computing Toolbox is not available; falling back to a serial loop.');
        cached = false;
    end
end
tf = cached;
end

% ----------------------------------------------------------------
function tf = canUseGPU()
persistent cached
if isempty(cached)
    try
        cached = gpuDeviceCount() > 0;
        if ~cached
            warning('backProjectionCore:noGPU', 'cfg.useGPU=true but no usable GPU was found; falling back to CPU.');
        end
    catch
        warning('backProjectionCore:noGPU', ...
            'cfg.useGPU=true but the Parallel Computing Toolbox / a supported GPU is not available; falling back to CPU.');
        cached = false;
    end
end
tf = cached;
end
