function cfg = config(mode)
% CONFIG  Central configuration for the Wide-FoV curved-array SISO-SAR
%          simulation (Wu, Gan & Chen, IEEE SPL 2025).
%
% cfg = config()          -> full paper-scale parameters (Section IV)
% cfg = config('quick')   -> small parameters for fast automated testing
%
% All numeric values below are taken directly from Section IV of the
% paper unless a comment says otherwise. Nothing else in the project
% should hard-code these numbers -- change them here.

if nargin < 1
    mode = 'full';
end
cfg.mode = mode;

cfg.c = 299792458;                 % speed of light, m/s

% ---------------------------------------------------------------------
% Proposed curved-linear-array system (Section IV, paragraph 1)
% ---------------------------------------------------------------------
cfg.curved.R0            = 0.3;    % curved array radius, m
cfg.curved.phi_c_deg      = 60;    % array central angle (max elevation), deg
cfg.curved.N              = 188;   % number of Tx/Rx pairs along the arc
cfg.curved.azimuth_step_deg = 0.38;% approx. azimuth sampling interval, deg
% M is not stated explicitly; the paper only gives the ~0.38 deg interval
% for 360 deg of rotation. We derive M from that interval (documented
% assumption -- see README "Limitations").
cfg.curved.M = round(360 / cfg.curved.azimuth_step_deg);

% ---------------------------------------------------------------------
% Planar-scanning comparison system (Section IV, paragraph 1)
% ---------------------------------------------------------------------
cfg.planar.length      = 0.6;      % linear array length, m
cfg.planar.numElements = 353;      % elements along the array
cfg.planar.numScan     = 353;      % mechanical scan positions

% ---------------------------------------------------------------------
% Frequency sweep (Section IV, paragraph 1)
% ---------------------------------------------------------------------
cfg.freq.f_start = 80e9;           % Hz
cfg.freq.f_end   = 96e9;           % Hz
cfg.freq.numFreq = 401;

% ---------------------------------------------------------------------
% Antenna
% ---------------------------------------------------------------------
cfg.HPBW_deg = 60;                 % half-power beamwidth, deg

% ---------------------------------------------------------------------
% Default single point-target location (used by experimentPointTarget)
% ---------------------------------------------------------------------
cfg.target.r         = 1.0;        % m
cfg.target.theta_deg  = 45;        % deg
cfg.target.phi_deg    = 30;        % deg

% Resolution-vs-elevation sweep (Fig. 5 methodology)
cfg.resSweep.r         = 1.0;
cfg.resSweep.theta_deg = 45;
cfg.resSweep.phi_deg_list = 0:5:90; % coarser than a continuous sweep for runtime

% ---------------------------------------------------------------------
% Imaging grid defaults (NOT given by the paper -- chosen for a
% tractable, honest demonstration; see README "Limitations")
% ---------------------------------------------------------------------
cfg.imaging.voxelRes         = 0.01;  % m, isotropic grid spacing
cfg.imaging.windowHalfWidth  = 0.15;  % m, half-width of local 3-D window
cfg.imaging.dB_range         = 15;    % dB, display floor (matches paper figures)

% ---------------------------------------------------------------------
% Reconstruction engine / performance controls
%   [EXISTING IMPLEMENTATION PARAMETERS -- not from the paper. These are
%    the original low-level controls and remain the authoritative knobs
%    used by backProjectionCore.m. The cfg.performance section further
%    down is a convenience layer that RESOLVES onto these; with the
%    default cfg.performance settings it leaves every value here
%    unchanged, so the original numerical behaviour is reproduced.]
% ---------------------------------------------------------------------
cfg.recon.method        = 'fft';   % 'fft' (fast, default) or 'direct' (brute force Eq.(8), validation only)
cfg.recon.zeroPadFactor = 8;       % range-profile oversampling for interpolation accuracy
cfg.recon.chunkSize     = 3000;    % aperture samples processed per chunk (memory bound)
cfg.useGPU              = false;   % optional; requires Parallel Computing Toolbox + supported GPU
cfg.useParallel         = false;   % optional; requires Parallel Computing Toolbox

% =====================================================================
% ADAPTIVE RECONSTRUCTION EXTENSIONS
%   [OUR PROPOSED EXTENSION PARAMETERS -- NOT from the paper. The paper
%    (Wu, Gan & Chen, IEEE SPL 2025) uses ONLY a fixed HPBW gate and the
%    geometry-only 1/count weight. Everything in this section is our own
%    contribution; see ADAPTIVE_METHODS.md for the equations, rationale
%    and limitations. All defaults below reproduce the ORIGINAL paper
%    behaviour exactly: 'original' weighting == 1/count, 'fixed' gating
%    == the nominal cfg.HPBW_deg gate.]
% =====================================================================

% ---- Adaptive spatial weighting (see computeAdaptiveWeight.m) --------
cfg.weighting.mode      = 'original';  % 'original' (Eq.9, 1/count) | 'adaptive'
cfg.weighting.alpha     = 0.5;         % confidence exponent; 0 -> reduces to 'original'
cfg.weighting.minFactor = 0.5;         % lower clip on the confidence factor (finite, >0)
cfg.weighting.maxFactor = 2.0;         % upper clip on the confidence factor (no runaway boost)
cfg.weighting.epsilon   = 1e-3;        % numerical floor for the normalised-magnitude ratio

% ---- Adaptive HPBW gating (see computeAdaptiveHPBW.m) ---------------
cfg.gating.mode          = 'fixed';    % 'fixed' (nominal cfg.HPBW_deg gate) | 'adaptive'
cfg.gating.minHPBW_deg   = 45;         % hard lower bound on the effective per-sample HPBW
cfg.gating.maxHPBW_deg   = 75;         % hard upper bound on the effective per-sample HPBW
cfg.gating.adaptStrength = 0.30;       % 0 -> reduces to 'fixed'; fraction of nominal HPBW
                                        % swung between central and peripheral rings

% ---------------------------------------------------------------------
% Noise / SNR (optional complex AWGN; see utils/addComplexAWGN.m)
%   [OUR EXTENSION -- the paper's Eq.(1)/(2) echo model is noiseless.
%    Used to show adaptive weighting is not merely amplifying noise.]
% ---------------------------------------------------------------------
cfg.noise.enable = false;   % default OFF -> echo generation is byte-identical to the original
cfg.noise.SNR_dB = 20;      % complex AWGN level (dB, relative to mean echo power)
cfg.noise.seed   = 0;       % RNG seed for reproducibility ([] -> do not reseed)

% ---------------------------------------------------------------------
% Runtime / performance execution framework
%   [OUR SECONDARY (implementation) CONTRIBUTION -- a convenience layer
%    over the existing cfg.recon.* / cfg.useGPU / cfg.useParallel knobs.
%    See PERFORMANCE_OPTIMIZATION.md. FFT evaluation, chunking, parfor
%    and gpuArray support ALREADY EXISTED in backProjectionCore.m; our
%    contribution is the unified mode selector, memory-aware automatic
%    chunk sizing, the auto hardware selection, the execution
%    diagnostics and the quantitative benchmark.]
%
%   mode: 'direct' | 'fft_serial' | 'fft_parallel' | 'fft_gpu' | 'fft_auto'
%   The DEFAULT 'fft_serial' resolves to method='fft', useGPU=false,
%   useParallel=false -- i.e. exactly the original defaults above.
% ---------------------------------------------------------------------
cfg.performance.mode           = 'fft_serial';
cfg.performance.chunkMode      = 'fixed';  % 'fixed' | 'auto'
cfg.performance.chunkSize      = 0;        % 0 -> keep cfg.recon.chunkSize; >0 -> override it
cfg.performance.targetMemoryMB = 512;      % memory budget per chunk for 'auto' chunk sizing
cfg.performance.verbose        = false;    % print requested/actual mode, #workers, chunk size

% ---------------------------------------------------------------------
% Quick-mode overrides: small, fast parameters for automated
% self-verification. Physically the same model, just far fewer samples.
% ---------------------------------------------------------------------
if strcmpi(mode, 'quick')
    % NOTE: N=24 / step=6deg (1440 raw samples) runs but is too sparse
    % for the artifact-suppression / uniformity benefit of the proposed
    % algorithm to show up in the point-array metrics (verified: it
    % actually got slightly worse at that density). N=48 / step=3deg
    % (5760 raw samples) reliably reproduces the paper's qualitative
    % trend while still running in a few seconds -- see
    % experiments/experimentPointArray.m and README "Limitations".
    cfg.curved.N = 48;
    cfg.curved.azimuth_step_deg = 3;
    cfg.curved.M = round(360 / cfg.curved.azimuth_step_deg);

    cfg.planar.numElements = 31;
    cfg.planar.numScan     = 31;

    % NOTE: numFreq cannot be reduced arbitrarily. The FFT-based
    % reconstruction (see backProjectionCore.m) has an unambiguous
    % range of c/(2*df), where df = bandwidth/(numFreq-1). Shrinking
    % numFreq while keeping the 80-96 GHz bandwidth fixed *increases*
    % df and *shrinks* the unambiguous range; if it falls below the
    % actual array-to-target distances (~0.7-2.2 m in these
    % experiments) the reconstruction aliases. 301 points keeps the
    % unambiguous range above 2.8 m while still being much faster to
    % IFFT than the full 401. validate_project.m checks this margin.
    cfg.freq.numFreq = 301;

    cfg.imaging.voxelRes        = 0.02;
    cfg.imaging.windowHalfWidth = 0.10;

    cfg.recon.zeroPadFactor = 4;
    cfg.recon.chunkSize     = 1500;

    cfg.resSweep.phi_deg_list = 0:15:90;

    % Our extensions: smaller per-chunk memory budget for the fast tests.
    cfg.performance.targetMemoryMB = 128;
end

% ---------------------------------------------------------------------
% Paths
% ---------------------------------------------------------------------
% Each call to config() starts a new timestamped run folder under
% results/, so results from different runs never overwrite each other.
% All experiments called with the SAME cfg (e.g. everything under one
% main.m session) share this one folder.
cfg.projectRoot = fileparts(mfilename('fullpath'));
cfg.runId       = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
cfg.resultsDir  = fullfile(cfg.projectRoot, 'results', ['run_' cfg.runId]);
if ~exist(cfg.resultsDir, 'dir')
    mkdir(cfg.resultsDir);
end

end
