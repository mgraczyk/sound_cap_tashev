%*******************************************
%
%   Selects the best microphone for given direction
%
%   Parameters:
%
%   MicArrDescriptor - matrix with the mcirophone array geometry
%   iAngle, Distance, iElev - desired listening direction
%                    angles in degrees, distance in meters
%
%   August 2005       Ivan Tashev
%
%*******************************************

function bestMic = BestMicrophone(MicArrDescriptor, iAngle, Distance, iElev)

if (nargin < 4)
    iElev = 0.0;
end

% constants
TORAD = pi / 180.0;
soundSpeed = 342;   % m/s

% indexes in MicArrDescriptor
cX = 1;
cY = 2;
cZ = 3;
aDir = 4;
aElev = 5;
cType = 6;
szMics = size(MicArrDescriptor);
nMics = szMics(1);

freq = 1000;

% Compute the gains
Angle = iAngle * TORAD;
sX = Distance * cos(Angle)*cos(iElev);
sY = Distance * sin(Angle)*cos(iElev);
sZ = Distance * sin(iElev);
    
for iMic = 1:nMics
    dist = sqrt((sX - MicArrDescriptor(iMic, cX))*(sX - MicArrDescriptor(iMic, cX)) + (sY - MicArrDescriptor(iMic, cY))*(sY - MicArrDescriptor(iMic, cY)) + (sZ - MicArrDescriptor(iMic, cZ))*(sZ - MicArrDescriptor(iMic, cZ)));
    gamma(iMic) = atan2((sY - MicArrDescriptor(iMic, cY)),(sX - MicArrDescriptor(iMic, cX)));
    gamma(iMic) = gamma(iMic) - MicArrDescriptor(iMic, aDir)*TORAD;
    distP(iMic) = sqrt((sX - MicArrDescriptor(iMic, cX))*(sX - MicArrDescriptor(iMic, cX)) + (sY - MicArrDescriptor(iMic, cY))*(sY - MicArrDescriptor(iMic, cY)));
    cappa(iMic) = atan2((sZ - MicArrDescriptor(iMic, cZ)), distP(iMic));
    cappa(iMic) = cappa(iMic) - MicArrDescriptor(iMic, aElev)*TORAD;
    CosTheta(iMic) = cos(gamma(iMic))*cos(cappa(iMic));
    Gain(iMic) = (1/dist) * exp(-j*2.0*pi*freq * dist / soundSpeed) * MicrophoneDirectivity(Distance, freq, CosTheta(iMic), MicArrDescriptor(iMic, cType));
end

[dMax, bestMic] = max(abs(Gain));

end
