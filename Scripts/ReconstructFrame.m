%==================================================
%
%   [Signal, prevFrameNext] = ReconstructFrame(Spec, prevFrame)
%
%   Reconstruct signal frame from the current frame spectra 
%               and previous frame
%
%   Spec          -   audio frame, frequency domain, length N
%   prevFrame     -   previous frame, time domain, length 2N
%   Signal        -   output audio frame, time domain, length N
%   prevFrameNext -   current frame, time domain, length 2N
%
%   (c) 2007 Ivan Tashev
%
%==================================================
function [Signal, prevFrameNext] = ReconstructFrame(Spec, prevFrame);

% restore the upper part of the frame spectrum
nSamples = length(Spec);
frameLen = 2*nSamples;
Spectr(1:nSamples) = Spec;
Spectr(nSamples+1) = 0.0 + i * 0.0;
Spectr(nSamples+2:2*nSamples) = conj(Spectr(nSamples:-1:2));

% compute the weighted time domain frame
t = (0:frameLen-1)/(frameLen-1);
ha = sin(pi*t);
prevFrameNext = ha .* ifft(Spectr);

% reconstruct using the second half of the previous frame  
% and the first half of the current frame
Signal = prevFrame(nSamples+1:2*nSamples)+prevFrameNext(1:nSamples);

