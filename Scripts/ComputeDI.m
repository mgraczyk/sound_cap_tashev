%==============================
%
%   Computes DI based on the Gain as function of the incident angle
%
% Parameters:
%
%   Angles   - array of angles from 0 to 180, in degrees
%   RMSgain  - in relative units, for the coresponding angle
%   angleMRA - main responce axis, in degrees
%
%   The program assumes simetry arround the main responce axis
%
%   (c) 2004  Ivan Tashev
%
%==============================

function DI = ComputeDI(Angles, RMSgain, angleMRA)

if (nargin < 3)
    angleMRA = 0;
end

DI = 0;
TORAD = pi/180;

numAngles = length(Angles);
numRMS = length(RMSgain);
if (numAngles ~= numRMS)
    disp('The vectors sizes should match!');
    return;
end
    
PowerIntegral = 0;
nPoints = 0;
Step = 5;

for fi = -180:Step:180
    for theta = 0:Step:90-Step
        dirAngle = NormAngle(fi - angleMRA);
        incidentAngle = sign(dirAngle) * acos(cos(dirAngle*TORAD)*cos(theta*TORAD))/TORAD;
        actualAngle = NormAngle(incidentAngle + angleMRA);
        RMS = interp1(Angles,RMSgain,actualAngle,'linear','extrap');
        PowerIntegral = PowerIntegral + RMS*RMS*cos(theta*TORAD);
        nPoints = nPoints + 1;
   end
end

PowerIntegral = PowerIntegral / nPoints;
RMS0 = interp1(Angles,RMSgain,angleMRA,'linear','extrap');

if (RMS0 > 0.0) && (PowerIntegral > 0.0)
    DI = RMS0*RMS0 / PowerIntegral;
    DI = 10 * log10(DI/1.50);
else
    DI = 0.0;
end

        
        
        