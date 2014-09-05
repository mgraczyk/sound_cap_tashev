%
%   Does C-weigting of given array for given set of frequencies
%
%   Arguments: Arrray - values for weigting
%              Freq   - frequency points
%
function WArray = Cweighting(Freq, Array)
S = j * 2.0 * pi * Freq;
GA = 5.8760e+009 * (S.^2) ./ ((S + 129.4336).^2) ./ ((S + 7.6655e+004).^2);
WArray = abs(GA) .* Array;

