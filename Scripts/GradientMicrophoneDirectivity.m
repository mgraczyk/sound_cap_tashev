%==================================================
%
%   MicGain = GradientMicrophoneDirectivity(Dist, Freq, CosTheta)
%
%   Generates the complex gain of a gradient microphone
%   given frequency, distance and incident angle cosine
%
%   Dist    -   distance from the sound source
%   Freq    -   signal frequency
%   CosTheta -  cosine of the incident angle
%
%   MicGain -   complex microphone gain
%
%   (c)2005 Ivan Tashev
% 
%==================================================

function MicGain = GradientMicrophoneDirectivity(Dist, Freq, CosTheta)

%
%   constants
%
SoundSpeed = 342;   % ms/ at 20 C and 1 MPa
Delta = 0.005;      % m microphone length
Tau = Delta / SoundSpeed;
omega = 2.0*pi*Freq;
omega1000 = 2.0*pi*1000;

%
%   ideal gradient microphone
%
if (Dist > 100.0*Delta)
    %   far field simple model
    cRatio = 1.0 - exp(-j * omega * Tau * CosTheta);
    cRatioMRA1000 = 1.0 - exp(-j * omega1000 * Tau);
else
    %   near field model
    rhoF = Delta / 2.0 + sqrt((Delta/2.0)^2+Dist^2-Dist*Delta*CosTheta);
    rhoB = Delta / 2.0 + sqrt((Delta/2.0)^2+Dist^2+Dist*Delta*CosTheta);
    cRatio = exp(-j*omega*rhoF/SoundSpeed) / rhoF - exp(-j*omega*rhoB/SoundSpeed) / rhoB;

    rhoFMRA = Dist;
    rhoBMRA = Dist + Delta;
    cRatioMRA1000 = exp(-j*omega1000*rhoFMRA/SoundSpeed) / rhoFMRA - exp(-j*omega1000*rhoBMRA/SoundSpeed) / rhoBMRA;
end

% normalize @ 1000 Hz and 0 degrees
cRatio = cRatio / cRatioMRA1000;

%
%   compensation filter 
%
%   filter paramters
FiltTau = 1e-3;
Cap = 33e-6;
Res = FiltTau / Cap;

% compute the capacitor impedance
if (Freq < 1e-10)
    Zc = - j * 1e-10;
else
    Zc = - j / (omega * Cap);
end
Zc1000 = - j / (omega1000 * Cap);

% compute the filter transfer gain
cGain = Zc / (Zc + Res);
cGain1000 = Zc1000 / (Zc1000 + Res);

%   normalize at 1000 Hz
cGain = cGain / cGain1000;

%
% comnbine
%
MicGain = cRatio * cGain;
