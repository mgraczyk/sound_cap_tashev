%==================================================
%
%   ambNoise = AmbientNoise(Freq)
%
%   Returns estimated magnitude of the ambient noise
%   given frequency, uses model from file NoiseModel.DAT.
%   The data file has two numbers on each line: frequency 
%   and amplitude. 
%
%   Freq    -   signal frequency
%
%   ambNoise -   real noise magnitude
%
%   (c)2005 Ivan Tashev
% 
%==================================================

function ambNoise = AmbientNoise(Freq)

persistent AmbientNoiseModel

if (isempty(AmbientNoiseModel))
    AmbientNoiseModel = load('NoiseModel.dat');
end

ambNoise = interp1(AmbientNoiseModel(:,1), AmbientNoiseModel(:,2), Freq, 'linear','extrap');

return;


