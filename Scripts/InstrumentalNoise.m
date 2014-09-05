%==================================================
%
%   instNoise = InstrumentalNoise(Freq)
%
%   Returns estimated magnitude of the instrumenta noise
%   given frequency, uses model from file AmpNoiseModel.DAT.
%   The data file has two numbers on each line: frequency 
%   and amplitude. 
%
%   Freq    -   signal frequency
%
%   instNoise -   real instrumental noise magnitude
%
%   (c)2005 Ivan Tashev
% 
%==================================================

function instNoise = InstrumentalNoise(Freq)

persistent InstrumentalNoiseModel;

if (isempty(InstrumentalNoiseModel))
    InstrumentalNoiseModel = load('AmpNoiseModel.dat');
end

instNoise = interp1(InstrumentalNoiseModel(:,1), InstrumentalNoiseModel(:,2), Freq, 'linear','extrap');

return;


