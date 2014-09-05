%==============================
%
% Weights file writing helper 
%
% Parameters:
%
%   fileName  - weights file name
%   samplRate - sampling rate in Hz
%   nMics     - number of microphones
%   frameSize - number of frequency bins
%   Weights   - complex matrix with the weights, 
%               size nBeams x nMics x numAngles
%   Angles    - angles the beams point to
%
%   (c) 2004  Ivan Tashev
%
%==============================

function WriteWeights(fileName, samplRate, nMics, frameSize, Weights, Angles)

WeightsRe = real(Weights);
WeightsIm = imag(Weights);
numBeams = length(Angles);

Elev = 0.0;
Dist = 1.5;
TORAD = pi / 180;

fid = fopen(fileName, 'w');
fprintf(fid,'%f\r\n', samplRate);       % sampling rate
fprintf(fid,'%d\r\n', nMics);           % # of channels
fprintf(fid,'%d\r\n', frameSize);       % # of coefficients
fprintf(fid,'%d\r\n\r\n', numBeams);    % # of beams

for iBeam = 1:numBeams
    fprintf(fid,' %16.8e %16.8e  %16.8e\r\n', Angles(iBeam)*TORAD, Elev, Dist);    % beam direction
    for iBin = 1:frameSize
        for iChan = 1:nMics
            fprintf(fid,'%18.8e %18.8e', WeightsRe(iBeam, iChan, iBin), WeightsIm(iBeam, iChan, iBin));
        end
        fprintf(fid,'\r\n');
    end
    fprintf(fid,'\r\n');
end
fclose(fid);

