function [ksdensity_return] = wrapped_ksdensity(ResidenceFrames, ResidenceFramesGoneFreq)
%% This is a wrapped ksdensity to wrap ksdensity return into one struct
% Zuhui Wang 2020/04/09
% Purpose: wrap up ksdensity returned value into one struct to be a valid
% function handle of bootstrp() to return both f and xi simultaneously

% Input:
%   ResidenceFrames: the filtered ResidenceFrames from track mat file
%   ResidenceFramesGoneFreq: the proper spaced x value (residence frames) of 1-CDF
% Output:
%   ksdensity_return.f: 1-CDF value at the corresponding ksdensity_return.xi
%   ksdensity_return.xi: see explanation above
%  ResidenceFramesGoneFreq{1} is pts for ksdensity 

[f, xi] = ksdensity(ResidenceFrames,ResidenceFramesGoneFreq{1},'Function','survivor');
ksdensity_return.f = f;
ksdensity_return.xi = xi;
end

