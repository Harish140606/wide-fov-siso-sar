function experimentSampling(cfg)
% EXPERIMENTSAMPLING  Part 14 / Fig. 4: visualize the raw (non-uniform)
% curved-array sampling aperture and the uniformized sub-aperture
% produced by Eq. (4)-(7).

if nargin < 1 || isempty(cfg)
    cfg = config();
end

fprintf('\n=== Experiment 1: Sampling Distribution ===\n');

ap = generateCurvedArray(cfg);          % Eq. (3)
subAp = createUniformSubAperture(ap);   % Eq. (4)-(7)

fprintf('Raw aperture samples (N=%d elevations x M=%d azimuths): %d\n', ap.N, ap.M, numel(ap.x));
fprintf('Uniform sub-aperture samples: %d (%.1f%% of raw)\n', ...
    numel(subAp.idx), 100*numel(subAp.idx)/numel(ap.x));
fprintf('Reference sampling distance d (Eq. 4): %.4f m\n', subAp.d);

plotSamplingPoints(ap.pos, ...
    sprintf('Raw Curved-Array Sampling (N=%d, M=%d, %d points)', ap.N, ap.M, numel(ap.x)), ...
    fullfile(cfg.resultsDir, 'sampling_before.png'));

plotSamplingPoints(subAp.pos, ...
    sprintf('Uniformized Sub-Aperture (%d points)', numel(subAp.idx)), ...
    fullfile(cfg.resultsDir, 'sampling_after.png'));

% Quantitative check: azimuthal sample count should now track sin(phi_i)
% roughly linearly (Eq. 5) rather than being constant M for every ring.
save(fullfile(cfg.resultsDir, 'sampling_data.mat'), 'ap', 'subAp');

fprintf('Saved: sampling_before.png, sampling_after.png, sampling_data.mat\n');
end
