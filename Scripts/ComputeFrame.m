%==================================================
%
%   Spectr = ComputeFrame(Signal)
%
%   Compute weighted short term spectrum of audio frame
%
%   Signal  -   audio frame, time domain, length 2N
%   Spectr  -   lower half of the spectrum, length N
%
%   (c) 2007 Ivan Tashev
%
%==================================================
function Spectr = ComputeFrame(Signal)

% weight and compute the spectrum
frameLen = length(Signal);
t = (0:frameLen-1)'/(frameLen-1);
ha = sin(pi*t);
Spec = fft(ha .* Signal);

% cut and return the lower half
Spectr = Spec(1:frameLen/2);

