%==================================================
%
%   MicGain = MicrophoneDirectivity(Dist, Freq, CosTheta, Type)
%
%   Generates the complex gain of a generic microphone
%   given frequency, distance and incident angle cosine
%
%   Dist    -   distance from the sound source
%   Freq    -   signal frequency
%   CosTheta -  cosine of the incident angle
%   Type    -   the microphone type, 
%               0 - omnidirectional (pressure)
%               1 - subcardiod
%               2 - cardiod
%               3 - spercardiod
%               4 - hypercardiod
%               5 - figure-8 (pressure gradient)
%
%   MicGain -   complex microphone gain
%
%   (c)2005 Ivan Tashev
% 
%==================================================
function MicGain = MicrophoneDirectivity(Dist, Freq, CosTheta, Type)

    switch Type
        case 0          % omnidirectional microphone
            Alpha = 1.0;
            Beta = 0.0;
        case 1          % subcardioid microphone
            Alpha = 0.75;
            Beta = 0.25;
        case 2          % cadrioid microphone
            Alpha = 0.5;
            Beta = 0.5;
        case 3          % supercardioid microphone
            Alpha = 0.37;
            Beta = 0.63;
        case 4          % hypercardioid microphone
            Alpha = 0.25;
            Beta = 0.75;
        case 5          % pressure gradient microphone
            Alpha = 0.0;
            Beta = 1.0;
    end

    MicGain = Alpha + Beta * GradientMicrophoneDirectivity(Dist, Freq, CosTheta);
    
end

