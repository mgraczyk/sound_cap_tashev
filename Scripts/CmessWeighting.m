%
%   Does CmessWeighting-weigting of given array for given set of frequencies
%
%   Arguments: Arrray - values for weigting
%              Freq   - frequency points
%
function WArray = CmessWeighting(Freq, Array)

FCmess = [    60    100    200    300    400   500   600   700   800   900 1000  1100  1200  1300  1400  1500  1600  1700  1800  1900  2000  2100  2200  2300  2400  2500  2600  2700  2800  2900  3000  3100  3200  3300  3400  3500  3600   3700   3800   3900   4000];
WCmess = [-54.65 -41.71 -25.17 -16.64 -11.29 -7.55 -4.75 -2.66 -1.19 -0.32 0.03  0.03 -0.17 -0.44 -0.71 -0.94 -1.12 -1.24 -1.32 -1.36 -1.38 -1.39 -1.41 -1.44 -1.50 -1.60 -1.76 -1.97 -2.26 -2.62 -3.09 -3.66 -4.35 -5.18 -6.18 -7.36 -8.75 -10.36 -12.12 -13.72 -14.43];

fMin = min(FCmess);
fMax = max(FCmess);
nFreq = length(Freq);
GA(1:nFreq) = 0.0;
for iFreq = 1:nFreq
    if (Freq(iFreq) >= fMin) && (Freq(iFreq) <= fMax)
        GA(iFreq) = 10.0^(interp1(FCmess,WCmess,Freq(iFreq))/20);
    end
end

WArray = abs(GA) .* Array;

