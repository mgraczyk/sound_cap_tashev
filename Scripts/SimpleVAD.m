%==============================
%
%   Simple voice activity detector
%
%   Parameters:
%
%   dLevel  - r.m.s. of the currecnt frame
%   dTime   - current time in seconds
%   VADdata - data sturcture with following memebrs:
%               tTauDown  - time constant, recomended value 40 ms
%               tTauUp    - time constant, recomended value 10 s
%               tThSpeech - threshold to speech, recomended value 3.5
%               tThNoise  - threshold to noise, recomended value 1.5
%               tLevel    - current minimal level, initialize with high value
%               tTime     - last time called, set 0 at the beginning
%               tSignal   - VAD state, 0 - noise frame, 1 - speech frame
%                           initialize with 0
%
%   (c) 2004  Ivan Tashev
%
%==============================
function VADdata = SimpleVAD(dLevel, dTime, VADdata)

if (nargin == 0)
    VADdata.tTauDown = 0.04;
    VADdata.tTauUp = 10.0;
    VADdata.tThSpeech = 3.5;
    VADdata.tTauDown = 0.04;
    VADdata.tThNoise = 1.5;
    VADdata.tLevel = 1000.0;
    VADdata.tTime = 0.0;
    VADdata.tSignal = 0;
    return;
end

if (dLevel < VADdata.tLevel)
    VADdata.tLevel = VADdata.tLevel + (dTime - VADdata.tTime) / VADdata.tTauDown * (dLevel - VADdata.tLevel);
    if (VADdata.tLevel < dLevel)
        VADdata.tLevel = dLevel;
    end
else
    VADdata.tLevel = VADdata.tLevel + (dTime - VADdata.tTime) / VADdata.tTauUp * (dLevel - VADdata.tLevel);
    if (VADdata.tLevel > dLevel)
        VADdata.tLevel = dLevel;
    end
end

if ((dLevel > VADdata.tThSpeech * VADdata.tLevel) && (dLevel > 1e-2))
    VADdata.tSignal = 1;
elseif (dLevel < VADdata.tThNoise * VADdata.tLevel)
    VADdata.tSignal = 0;
end

VADdata.tTime = dTime;
