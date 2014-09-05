%==============================
%
% Weights file reading helper 
%
% Parameters:
%
%   fileName - file name of the weights file
%
% Return paramteres:
%
%   samplRate - sampling rate in Hz
%   nMics     - number of microphones
%   frameSize - number of frequency bins
%   nBeams    - number of pre-computed beams
%   Weights   - complex matrix with the weights, 
%               size nBeams x nMics x numAngles
%   Angles    - angles the beams point to
%
%   (c) 2004  Ivan Tashev
%
%==============================

function [samplRate, nMics, frameSize, nBeams, Weights, Angles] = ReadWeights(fileName)

TORAD = pi / 180;

% Read the weights
fid = fopen(fileName, 'r');

strLine = fgets(fid);
samplRate = sscanf(strLine,'%f'); % sampling rate
strLine = fgets(fid);
nMics = sscanf(strLine,'%d');   % # of channels
strLine = fgets(fid);
frameSize = sscanf(strLine,'%d'); % # of coefficients
strLine = fgets(fid);
nBeams = sscanf(strLine,'%d');    % # number of beams
strLine = fgets(fid);

for iBeam = 1:nBeams
    strLine = fgets(fid);
    beamCoord = sscanf(strLine,'%f', [3]); % beam direction
    Angles(iBeam) = beamCoord(1) / TORAD;

    for iBin = 1:frameSize
        strLine = fgets(fid);
        WeightsRaw(1:2*nMics) = sscanf(strLine,'%e', [2*nMics]);
        for iMic = 1:nMics
            Weights(iBeam, iMic, iBin) = WeightsRaw(2*iMic-1) + i*WeightsRaw(2*iMic);
        end
    end
    strLine = fgets(fid);
end

fclose(fid);
