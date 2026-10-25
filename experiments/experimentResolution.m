function results = experimentResolution(cfg)
% EXPERIMENTRESOLUTION  Part 16 / Fig. 5 methodology: for a single-voxel
% scatterer at (r=1m, theta=45deg, phi in [0,90]deg), reconstruct the
% IPSF and measure the -3 dB widths of line slices along the azimuthal
% and range directions, for three system/algorithm combinations:
%   1) Planar-scanning SISO-SAR + classical BPA
%   2) Proposed curved SISO-SAR + classical BPA
%   3) Proposed curved SISO-SAR + proposed measurement-fusion algorithm
%
% Following the paper's own wording ("spatial resolutions are
% quantified by the -3 dB widths of IPSF LINE SLICES"), each resolution
% value is measured from a fine 1-D line of voxels through the target
% along the relevant direction, not from a full 3-D/2-D grid -- this
% matches the described methodology exactly and keeps the sweep fast
% enough to run automatically.

if nargin < 1 || isempty(cfg)
    cfg = config();
end

fprintf('\n=== Experiment 3: Resolution vs Elevation Angle ===\n');

curvedAp = generateCurvedArray(cfg);
planarAp = generatePlanarArray(cfg);
[~, k] = generateFrequencyGrid(cfg);

phiList = cfg.resSweep.phi_deg_list;
r = cfg.resSweep.r;
theta_deg = cfg.resSweep.theta_deg;

lineHalfWidth = 0.06;   % m, local window around the target for each 1-D cut
lineRes       = 0.001;  % m, fine sampling for accurate -3 dB crossing

nA = numel(phiList);
azRes  = nan(nA, 3);   % columns: [planar+BPA, curved+BPA, curved+proposed]
rgRes  = nan(nA, 3);
peakAmp = nan(nA, 3);  % ADDITIONAL diagnostic (not a paper-reported quantity, see README):
                        % raw (non peak-normalized) reconstructed amplitude. The antenna
                        % directivity taper (antennaGain.m) attenuates SNR/visibility at
                        % large off-boresight angles, but a peak-normalized -3dB width is,
                        % by construction, blind to that in a noiseless simulation -- this
                        % tracks the effect the -3dB width metric cannot see.

for a = 1:nA
    phi_deg = phiList(a);
    tgt = createPointTarget(r, theta_deg, phi_deg);

    theta_rad = deg2rad(theta_deg);
    rHat  = tgt.pos / norm(tgt.pos);
    azHat = [-sin(theta_rad), cos(theta_rad), 0];   % tangential (cross-range) direction

    tLine = (-lineHalfWidth:lineRes:lineHalfWidth)';
    rangeLine = tgt.pos + tLine * rHat;
    azLine    = tgt.pos + tLine * azHat;

    % ---- System 1: planar SISO-SAR + classical BPA ----
    echoP = generatePointTargetEcho(planarAp.pos, tgt.pos, tgt.amp, k, planarAp.beamDir, cfg.HPBW_deg);
    imgP_rg = BPA(echoP, planarAp.pos, k, rangeLine, cfg);
    imgP_az = BPA(echoP, planarAp.pos, k, azLine, cfg);
    rgRes(a,1) = computeResolution(tLine, abs(imgP_rg));
    azRes(a,1) = computeResolution(tLine, abs(imgP_az));
    peakAmp(a,1) = max(abs(imgP_az));

    % ---- System 2: curved SISO-SAR + classical BPA ----
    echoC = generatePointTargetEcho(curvedAp.pos, tgt.pos, tgt.amp, k, curvedAp.beamDir, cfg.HPBW_deg);
    imgC_rg = BPA(echoC, curvedAp.pos, k, rangeLine, cfg);
    imgC_az = BPA(echoC, curvedAp.pos, k, azLine, cfg);
    rgRes(a,2) = computeResolution(tLine, abs(imgC_rg));
    azRes(a,2) = computeResolution(tLine, abs(imgC_az));
    peakAmp(a,2) = max(abs(imgC_az));

    % ---- System 3: curved SISO-SAR + proposed measurement fusion ----
    imgM_rg = measurementFusionReconstruction(curvedAp, echoC, k, rangeLine, cfg);
    imgM_az = measurementFusionReconstruction(curvedAp, echoC, k, azLine, cfg);
    rgRes(a,3) = computeResolution(tLine, abs(imgM_rg));
    azRes(a,3) = computeResolution(tLine, abs(imgM_az));
    peakAmp(a,3) = max(abs(imgM_az));

    fprintf('  phi=%5.1f deg | az[cm]: planar=%s curved-BPA=%s proposed=%s | rg[cm]: planar=%s curved-BPA=%s proposed=%s\n', ...
        phi_deg, fmt(azRes(a,1)), fmt(azRes(a,2)), fmt(azRes(a,3)), ...
        fmt(rgRes(a,1)), fmt(rgRes(a,2)), fmt(rgRes(a,3)));
end

legendCell = {'Planar SISO-SAR + BPA', 'Curved SISO-SAR + BPA', 'Curved SISO-SAR + Proposed'};
plotResolution(phiList, {azRes(:,1), azRes(:,2), azRes(:,3)}, legendCell, ...
    'Azimuthal Resolution vs Elevation Angle', 'resolution / cm', ...
    fullfile(cfg.resultsDir, 'resolution_azimuth.png'));
plotResolution(phiList, {rgRes(:,1), rgRes(:,2), rgRes(:,3)}, legendCell, ...
    'Range Resolution vs Elevation Angle', 'resolution / cm', ...
    fullfile(cfg.resultsDir, 'resolution_range.png'));

% ADDITIONAL diagnostic (not a paper-reported plot): raw peak amplitude
% vs angle, dB relative to each system's own phi=0 value. This is the
% quantity that actually captures the "planar system fails at large
% angles" effect in our noiseless simulation (see comment above).
peakAmpDB = 20*log10(peakAmp ./ (peakAmp(1,:) + eps) + eps);
plotResolution(phiList, {peakAmpDB(:,1)/100, peakAmpDB(:,2)/100, peakAmpDB(:,3)/100}, legendCell, ...
    'Reconstructed Peak Amplitude vs Elevation Angle (relative dB, ADDITIONAL diagnostic)', ...
    'relative amplitude / dB (not cm)', ...
    fullfile(cfg.resultsDir, 'resolution_peakAmplitude.png'));

results.phiList = phiList;
results.azRes = azRes;
results.rgRes = rgRes;
results.peakAmp = peakAmp;
results.legend = legendCell;
save(fullfile(cfg.resultsDir, 'resolution_data.mat'), 'results');

fprintf('Saved: resolution_azimuth.png, resolution_range.png, resolution_peakAmplitude.png, resolution_data.mat\n');
end

function str = fmt(v)
if isnan(v)
    str = 'NaN';
else
    str = sprintf('%.3f', v*100);
end
end
