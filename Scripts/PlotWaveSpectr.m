%==================================================
%
%   PlotWavSpectr(<inp_file_1> [,<inp_file_2>])
%
%   <inp_file_1>  -   audio file name in WAV format
%   <inp_file_2>  -   audio file name in WAV format, optional
%
%   (c) 2007 Ivan Tashev
%
%==================================================
function PlotWaveSpectr(fileName1, fileName2)

nFreqs = 512;

%
%   Plot the first signal
%
[signal1,fs1] = wavread(fileName1);
[Spectr1 Freqs1] = ComputeSpectra(signal1, fs1, nFreqs);
plot(Freqs1,Spectr1, 'b','LineWidth',2);
hold on

%
%   Plot the second signal if specified
%
if (nargin > 1)
    [signal2,fs2] = wavread(fileName2);
    [Spectr2 Freqs2] = ComputeSpectra(signal2, fs2, nFreqs);
    plot(Freqs2,Spectr2,'g','LineWidth',2);
    hold on
    legend(fileName1,fileName2);
else
    legend(fileName1);
end

%
%   Finalize the chart
%
hold off
title('Magnitude spectra','FontSize',14);
xlabel('Frequency, Hz','FontSize',14);
ylabel('Magnitude, dB','FontSize',14);
box('on');
grid('on')

return;
end

%
% compute average signal spectrum
%
function [Spectr Freq] = ComputeSpectra(Samples, Fs, nFreqs)

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
frameStep = nFreqs;
weightwindow = hann(frameSize); 

Spectr(1:nFreqs) = 0.0;
nFrames = 0;

for sample = 1:frameStep:nPoints-frameSize
    MomSpectr = abs(fft(weightwindow .* Samples(sample:sample+frameSize-1)));
    Spectr = Spectr + MomSpectr(1:nFreqs)';
    nFrames = nFrames + 1;
end

Spectr = 20.0 * log10(Spectr/nFrames);
return;
end

