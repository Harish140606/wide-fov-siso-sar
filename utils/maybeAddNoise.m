function echo = maybeAddNoise(echo, cfg)
% MAYBEADDNOISE  OUR EXTENSION. Convenience guard: applies complex AWGN
% to an echo matrix ONLY when cfg.noise.enable is true. With the default
% config (cfg.noise.enable = false) this returns the echo completely
% unchanged and does not even touch the RNG, so every existing
% experiment reproduces its original numerical result.
%
%   echo : Ns x Nf complex
%   cfg  : configuration struct (uses cfg.noise.enable / SNR_dB / seed)

if ~isfield(cfg, 'noise') || ~isfield(cfg.noise, 'enable') || ~cfg.noise.enable
    return;
end

seed = [];
if isfield(cfg.noise, 'seed')
    seed = cfg.noise.seed;
end
echo = addComplexAWGN(echo, cfg.noise.SNR_dB, seed);
end
