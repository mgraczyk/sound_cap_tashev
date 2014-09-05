function dRMS = ComputeRMS(Array, begBin, endBin)

if (nargin < 3)
    endBin = int32(length(Array));
end
if (nargin < 2)
    begBin = int32(1);
end

dRMS = sqrt(mean(abs(Array(begBin:endBin)).*abs(Array(begBin:endBin))));

%Arr = abs(Array);
%ArrSqr = Arr .* Arr;
%ArrAv = mean(ArrSqr);
%dRMS = sqrt(ArrAv);