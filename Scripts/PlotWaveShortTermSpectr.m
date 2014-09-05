%==================================================
%
%   PlotWavSpectr(fileName [,plotType])
%
%   Plots the spectrogram of a WAV file.
%
%   fileName  -   audio file name in WAV format
%   plotType  -   'S' for planar plot, 'D' for 3D plot
%
%   (c) 2007 Ivan Tashev
%
%==================================================
function PlotWaveShortTermSpectr(fileName, plotType)

if (nargin < 2)
    plotType = 'S';
end
plotType = upper(plotType);

% set some defaults
nFreqs = 512;

%
%   Compute the shot time magnitude spectrum
%
[signal,fs] = wavread(fileName);
[Spectr Freq Time] = ComputeShortTermSpectrum(signal, fs, nFreqs);

%
%   Plot the signal
%
minMag = -80.0;
maxMag = 0.0;

if (plotType == 'S')
    imagesc(Time,Freq,Spectr',[minMag maxMag]);
    axis xy;
    colorbar();
end
if (plotType == 'D')
    mesh(Time, Freq, min(maxMag,max(minMag,Spectr')));
    view([340,60])
    axis([min(Time) max(Time) min(Freq) max(Freq) -80 0]);
    zlabel('Magnitude, dB','FontSize',14);
end

xlabel('Time, sec','FontSize',14);
ylabel('Frequency, Hz','FontSize',14);
title(sprintf('Short-Time Magnitude Spectrum: %s',fileName),'FontSize',14);

return;

%
% compute short-time signal spectrum
%
function [Spectr Freq Time] = ComputeShortTermSpectrum(Samples, Fs, nFreqs)

if (nargin < 3)
    nFreqs = 512;
end

nPoints = length(Samples);
if (nPoints < 2*nFreqs)
    nFreqs = int32(nPoints / 2);
end
    
Freq = 1:nFreqs;
Freq = Freq * Fs / 2 / nFreqs;

frameSize = 2 * nFreqs;
frameStep = nFreqs/4;
weightwindow = hann(frameSize); 

nFrames = 1;
for sample = 1:frameStep:nPoints-frameSize
    MomSpectr = abs(fft(weightwindow .* Samples(sample:sample+frameSize-1)));
    Spectr(nFrames,:) = MomSpectr(1:nFreqs);
    nFrames = nFrames + 1;
end

Spectr = 20.0 * log10(max(1e-6,Spectr/nFrames));
Time = 1:nFrames-1;
Time = (Time - 0.5) * frameStep / Fs;
return;
