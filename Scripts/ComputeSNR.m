function dSNR = ComputeSNR(dSignalAndNoise, dNoise)

if ((dSignalAndNoise > dNoise) && (dNoise > 0))
    dSignal = sqrt(dSignalAndNoise.^2 - dNoise.^2);
    dS = 20 * log10(dSignal);
    dN = 20 * log10(dNoise);
    dSNR = dS - dN;
else
    dSNR = -80;
end

    
