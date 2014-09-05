%==================================================
%
%   [Signal, prevFrameNext] = RestoreFFT(Spec, prevFrame)
%
%   Reconstruct signal frame from the current frame spectra 
%               and previous frame
%
%   Spec          -   audio frame, frequency domain, length N
%   prevFrame     -   previous frame, time domain, length 2N
%   Signal        -   output audio frame, time domain, length N
%   prevFrameNext -   full frame, time domain, length 2N
%
%   (c) 2007 Ivan Tashev
%
%==================================================
function [Signal, prevFrameNext] = RestoreFFT(Spec, prevFrame)
    nSamples = length(Spec);
    Spectr(1:nSamples) = Spec;
    Spectr(nSamples+1) = 0.0 + i * 0.0;
    Spectr(nSamples+2:2*nSamples) = conj(Spectr(nSamples:-1:2));
    ha = sqrt(HANN(2*nSamples)');
    prevFrameNext = ha .* ifft(Spectr);
    Signal = prevFrame(nSamples+1:2*nSamples)+prevFrameNext(1:nSamples);
end

