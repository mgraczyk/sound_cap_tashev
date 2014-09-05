%*******************************************
%
%   Computes the noise gain of a microphone array
%
%   Parameters:
%
%   Weights          - weights matrix
%   MicArrDescriptor - matrix with the mcirophone array geometry
%   Freq             - frequency, Hz
%
%   August 2005       Ivan Tashev
%
%*******************************************

function NoiseGain = ComputeNoiseGain(Weights, MicArrDescriptor, Freq)

cX = 1;
cY = 2;
cZ = 3;
aDir = 4;
aElev = 5;
cType = 6;
szMics = size(MicArrDescriptor);
nMics = szMics(1);

TORAD = pi / 180.0;
soundSpeed = 342;   % m/s at 20 degrees Celsius, 1 MPa
Distance = 1.5; % m
Elev = 0.0; % elevation, degrees

NoiseGain = 0.0;
nGains = 0;

for dir = 0:10:350
    sX = Distance * cos(dir * TORAD) * cos(Elev * TORAD);
    sY = Distance * sin(dir * TORAD) * cos(Elev * TORAD);
    sZ = Distance * sin(Elev * TORAD);

    for iMic = 1:nMics
        dist(iMic) = sqrt((sX - MicArrDescriptor(iMic, cX))*(sX - MicArrDescriptor(iMic, cX)) + (sY - MicArrDescriptor(iMic, cY))*(sY - MicArrDescriptor(iMic, cY)) + (sZ - MicArrDescriptor(iMic, cZ))*(sZ - MicArrDescriptor(iMic, cZ)));
        gamma(iMic) = atan2((sY - MicArrDescriptor(iMic, cY)),(sX - MicArrDescriptor(iMic, cX)));
        gamma(iMic) = gamma(iMic) - MicArrDescriptor(iMic, aDir)*TORAD;
        distP(iMic) = sqrt((sX - MicArrDescriptor(iMic, cX))*(sX - MicArrDescriptor(iMic, cX)) + (sY - MicArrDescriptor(iMic, cY))*(sY - MicArrDescriptor(iMic, cY)));
        cappa(iMic) = atan2((sZ - MicArrDescriptor(iMic, cZ)), distP(iMic));
        cappa(iMic) = cappa(iMic) - MicArrDescriptor(iMic, aElev)*TORAD;
        CosTheta(iMic) = cos(gamma(iMic))*cos(cappa(iMic));
    end
    
    OmniGain = (1.0/Distance) * exp(-j*2.0*pi*Freq * Distance / soundSpeed);

    for iMic = 1:nMics
        Gain(iMic) = (1/dist(iMic)) * exp(-j*2.0*pi*Freq * dist(iMic) / soundSpeed);
        MicDir(iMic) = MicrophoneDirectivity(dist(iMic), Freq, CosTheta(iMic), MicArrDescriptor(iMic, cType));
        InpSignal(iMic) = Gain(iMic) * MicDir(iMic);
    end
        
    OutSignal = sum(InpSignal .* Weights);
    
    NoiseGain = NoiseGain + (abs(OutSignal) / abs(OmniGain))^2;
    nGains = nGains + 1;
end

AmbNoiseGain = sqrt(NoiseGain / nGains);

AmbNoise = AmbientNoise(Freq);
NCNoiseGain = sqrt(sum(abs(Weights).^2));
NCNoise = InstrumentalNoise(Freq);
NoiseGain = sqrt(((AmbNoise * AmbNoiseGain)^2 + (NCNoise * NCNoiseGain)^2) / (AmbNoise^2 + NCNoise^2));

end

