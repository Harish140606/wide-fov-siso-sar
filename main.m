% MAIN  Entry point for the Wide-FoV Curved SISO-SAR simulation project
% (numerical reproduction of Wu, Gan & Chen, "Wide Field-of-View MMW
% SISO-SAR Image Reconstruction Based on Curved Linear Array," IEEE
% Signal Processing Letters, vol. 32, 2025).
%
% Run this script (type `main` at the MATLAB prompt with this folder
% and its subfolders on the path) to get an interactive menu of
% experiments. See QUICK_START.md / USER_GUIDE.md for details.

if isempty(which('generateCurvedArray'))
    addpath(genpath(fileparts(mfilename('fullpath'))));
end

names = { ...
    'Sampling Distribution (Fig. 4)'; ...
    'Single Point Target (3-D IPSF)'; ...
    'Resolution vs Elevation Angle (Fig. 5)'; ...
    'Planar Point Array (Fig. 6)'; ...
    'Human Point Cloud (Fig. 8)'; ...
    'Adaptive Methods Comparison (our extension)'; ...
    'Runtime Benchmark (our extension)'; ...
    'Noise / SNR Robustness (our extension)'; ...
    'Run All Experiments'; ...
    'Run Validation Suite Only'; ...
    'Exit'};

fprintf('\n');
fprintf('===========================================\n');
fprintf(' Wide-FoV Curved SISO-SAR Simulation\n');
fprintf(' (Wu, Gan & Chen, IEEE SPL 2025 -- numerical reproduction)\n');
fprintf('===========================================\n');
for i = 1:numel(names)
    fprintf(' %d. %s\n', i, names{i});
end
fprintf('===========================================\n');
fprintf(' Configuration: type "cfg = config();" for full paper-scale\n');
fprintf(' parameters, or "cfg = config(''quick'');" for a fast demo.\n');
fprintf('===========================================\n');

choice = [];
while isempty(choice)
    raw = strtrim(input('Select an option [1-11] (you can also type the option name): ', 's'));
    if isempty(raw)
        choice = numel(names);   % default: Exit
        break;
    end
    numChoice = str2double(raw);
    if ~isnan(numChoice) && numChoice >= 1 && numChoice <= numel(names)
        choice = round(numChoice);
    else
        match = find(contains(lower(names), lower(raw)), 1);
        if ~isempty(match)
            choice = match;
        else
            fprintf('Not recognized: "%s". Enter a number 1-8 or part of an option name.\n', raw);
        end
    end
end
fprintf('-> Selected: %d. %s\n', choice, names{choice});

if choice == numel(names)
    fprintf('Exiting.\n');
    return;
end

modeChoice = input('Use (q)uick demo parameters or (f)ull paper-scale parameters? [q/f]: ', 's');
if isempty(modeChoice) || lower(modeChoice(1)) == 'q'
    cfg = config('quick');
    fprintf('Using QUICK configuration.\n');
else
    cfg = config('full');
    fprintf('Using FULL paper-scale configuration -- this can take a long time, see USER_GUIDE.md.\n');
end

switch choice
    case 1
        experimentSampling(cfg);
    case 2
        experimentPointTarget(cfg);
    case 3
        experimentResolution(cfg);
    case 4
        experimentPointArray(cfg);
    case 5
        experimentHuman(cfg);
    case 6
        experimentAdaptiveMethods(cfg);
    case 7
        experimentRuntimeBenchmark(cfg);
    case 8
        experimentNoiseRobustness(cfg);
    case 9
        experimentSampling(cfg);
        experimentPointTarget(cfg);
        experimentResolution(cfg);
        experimentPointArray(cfg);
        experimentHuman(cfg);
        experimentAdaptiveMethods(cfg);
        experimentRuntimeBenchmark(cfg);
        experimentNoiseRobustness(cfg);
    case 10
        validate_project();
end

fprintf('\nDone. Results saved under: %s\n', cfg.resultsDir);
