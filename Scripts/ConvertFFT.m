%==================================================
%
%   Spectr = ConvertFFT(Signal)
%
%   Compute weighted short term spectrum of audio frame
%
%   Signal  -   audio frame, time domain, length 2N
%   Spectr  -   lower half of the spectrum, length N
%
%   (c) 2007 Ivan Tashev
%
%==================================================
function Spectr = ConvertFFT(Signal)

frameLen = length(Signal);
ha = sqrt(HANN(frameLen));
Spec = fft(ha .* Signal);
Spectr = Spec(1:frameLen/2);

