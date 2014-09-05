%*******************************************
%
%   Microphone array analysing program
%
%   Parameters:
%
%   direction   - selects the closest beam
%                 defaults to 0
%   fileWeights - weights file name
%                 defaults to Weights.DAT
%
%   August 2005       Ivan Tashev
%
%*******************************************

function MAAnalysis(direction, fileWeights)

if (nargin < 1)
    direction = 0;
end
direction = NormAngle(direction);

if (nargin < 2)
    fileWeights = 'Weights.DAT';
end

%=============================================
%
%   Read the input data 
%
%=============================================

disp(' ');
disp('  Microphone array analysis');
disp('  (c) Ivan Tashev 2005');
disp('  =========================');
disp(' ');
disp(sprintf('Processing: %s', fileWeights));
disp(sprintf('Direction: %d degrees', direction));

% constants
TORAD = pi / 180.0;
soundSpeed = 342;   % m/s
Distance = 1.5; % m
Elev = 0.0; % elevation, degrees

MinGain = -25.0;
MaxGain = 3.0;

angleStep = 5;  % degrees
numAngles = 360 / angleStep + 1;
probeAngles = 0:angleStep:360;
probeAngles = probeAngles -180;
[dMin, indBeam] = min(abs(probeAngles-direction));

% Read the weights
[samplRate, nMics, frameSize, nBeams, Weights, Angles] = ReadWeights(fileWeights);
[dMin, curBeam] = min(abs(Angles - direction));

% Read the microphone array geometry
MicArrDescriptor = load('MicArrDesc.dat');
cX = 1;
cY = 2;
cZ = 3;
aDir = 4;
aElev = 5;
cType = 6;
szMics = size(MicArrDescriptor);
nDMics = szMics(1);
if (nMics ~= nDMics)
    disp('Different number of microphones!');
    return;
end
bestMic = BestMicrophone(MicArrDescriptor, direction, Distance);

disp('Input files reading completed');

%=============================================
%
%   Estimation loop
%
%=============================================

% Frequency response parameteres
freqStep = samplRate / frameSize / 2;
begFreq = freqStep / 2;
Frequency(1:frameSize) = begFreq:freqStep:samplRate/2;

ArrayGain(1:frameSize, 1:numAngles) = 1.0 + i * 0.0;

for iBin = 1:frameSize
    
    Freq = Frequency(iBin);
    
    for iAngle = 1:numAngles
        
        dir = probeAngles(iAngle);
        ArrayGain(iBin, iAngle) = ComputeArrayGain(dir, Elev, Distance, Freq, MicArrDescriptor, squeeze(Weights(curBeam, :, iBin)));
        
        MagnitudeLin(iBin, iAngle) = abs(ArrayGain(iBin, iAngle));
        if (MagnitudeLin(iBin, iAngle) > 0)
            Magnitude(iBin, iAngle) = 20.0 * log10(MagnitudeLin(iBin, iAngle));
        else
            Magnitude(iBin, iAngle) = -100.0;
        end
    end
end

Magnitude(1:frameSize, 1:numAngles) = max(Magnitude(1:frameSize, 1:numAngles), MinGain);
Magnitude(1:frameSize, 1:numAngles) = min(Magnitude(1:frameSize, 1:numAngles), MaxGain);

disp('Computing the response completed');

close all
fig = 0;

%
%   Plot the response - magnitude
%
fig = fig + 1;
[X, Y] = meshgrid(probeAngles, Frequency);
fig = fig + 1;
h = figure(fig);
mesh(X, Y, Magnitude);
view([240,40])
axis([-180 180 0 8000 MinGain MaxGain]);
xlabel('Angle, deg','FontSize',12);
ylabel('Frequency, Hz','FontSize',12);
zlabel('Gain, dB','FontSize',12);
title('Microphone Array Magnitude Response','FontSize',14);

disp('Plotting magnitude response completed');

%
%   Plot the response - phase
%
Phase = angle(ArrayGain) / TORAD;
fig = fig + 1;
h = figure(fig);
mesh(X, Y, Phase);
view([240,40])
axis([-180 180 0 8000 -180 180.0]);
xlabel('Angle, deg','FontSize',12);
ylabel('Frequency, Hz','FontSize',12);
zlabel('Phase, deg','FontSize',12);
title('Microphone Array Phase Response','FontSize',14);

disp('Plotting phase response completed');

%
%   Plot the response at various frequencies
%
AnglesRad = probeAngles * TORAD;
Freqs = [500, 1000, 2500, 5000];

for iFreq = 1: length(Freqs)
    freq = Freqs(iFreq);
    [val freqIndex] = min(abs(Frequency-freq));
    Response = squeeze(MagnitudeLin(freqIndex,:));
       
    fig = fig+1;
    h = figure(fig);
    polar(AnglesRad, Response);
    title(sprintf('Directivity pattern %5.0f Hz',freq),'FontSize',14);

    ResponseDb = squeeze(Magnitude(freqIndex,:));
    fig = fig + 1;
    h = figure(fig);
    plot(probeAngles, ResponseDb,'-b','LineWidth',2);
    axis([-180, 180, MinGain, max([max(Response), MaxGain])]);
    xlabel('Incident Angle, deg','FontSize',12);
    ylabel('Gain, dB','FontSize',12);
    title(sprintf('Directivity pattern %5.0f Hz',freq),'FontSize',14);
    grid('on')
end

%
%   Plot the response towards MRA
%
[dMin, indBeam] = min(abs(probeAngles - direction));
MagnitudeMRA = squeeze(Magnitude(:, indBeam));

fig = fig + 1;
h = figure(fig);
plot(Frequency, MagnitudeMRA,'-b','LineWidth',2);
axis([0, samplRate/2.0, MinGain, max([max(MagnitudeMRA), MaxGain])]);
xlabel('Frequency, Hz','FontSize',12);
ylabel('Gain, dB','FontSize',12);
title('Frequency Response','FontSize',14);
grid('on')

disp('Plotting MRA magnitude response completed');

%
%   Plot the directivity index
%
% compute DI per bin
for iBin = 1:frameSize
    DI(iBin) = ComputeDI(probeAngles, squeeze(MagnitudeLin(iBin, :)), direction);
end

% plot DI
fig = fig + 1;
h = figure(fig);
plot(Frequency,DI,'-b','LineWidth',2);
xlabel('Frequency, Hz','FontSize',12);
ylabel('DI, dB','FontSize',12);
title('Microphone Array Directivity Index','FontSize',14);
grid('on');

% print total DI
TotalDI = mean(DI);
disp(sprintf('Total Directivity Index: %6.2f dB', TotalDI));
DIlin = power(10, -DI/20);

%
%   Plot the noise suppression
%
for iBin = 1:frameSize
    Freq = Frequency(iBin);
    AmbNoise(iBin) = AmbientNoise(Freq);
    InstNoise(iBin) = InstrumentalNoise(Freq);
    OmniNoise(iBin) = sqrt(AmbNoise(iBin)*AmbNoise(iBin) + InstNoise(iBin)*InstNoise(iBin));
end

for iBin = 1:frameSize
    % compute the instrumental noise gain    
    InstrGain(iBin) = 0;
    for iMic = 1:nMics
        InstrGain(iBin) = InstrGain(iBin) + abs(Weights(curBeam, iMic, iBin)).^2;
    end
    InstrGain(iBin) = sqrt(InstrGain(iBin));

    % compute the ambient noise gain
    AmbientGain(iBin) = DIlin(iBin) * MagnitudeLin(iBin, indBeam);

    % compute the overal noise gain
    MANoise(iBin) = sqrt((AmbientGain(iBin)*AmbNoise(iBin))^2 + (InstrGain(iBin)*InstNoise(iBin))^2);
    
    if (OmniNoise(iBin) > 0.0)
        NoiseGain(iBin) = MANoise(iBin) / (OmniNoise(iBin));
    else
        NoiseGain(iBin) = 0.0;
    end
    
    if (NoiseGain(iBin) > 0.00001)
        NoiseGainL(iBin) = 20.0 * log10(NoiseGain(iBin));
    else
        NoiseGainL(iBin) = -100.0;
    end
    
    if (AmbientGain(iBin) > 0.00001)
        AmbientGainL(iBin) = 20.0 * log10(AmbientGain(iBin));
    else
        AmbientGainL(iBin) = -100.0;
    end
    
    if (InstrGain(iBin) > 0.00001)
        InstrGainL(iBin) = 20.0 * log10(InstrGain(iBin));
    else
        InstrGainL(iBin) = -100.0;
    end
        
end

fig = fig + 1;
h = figure(fig);
plot(Frequency,NoiseGainL,'-b','LineWidth',2);
hold('on');
plot(Frequency,AmbientGainL,'-r','LineWidth',2);
hold('on');
plot(Frequency,InstrGainL,'-g','LineWidth',2);
hold('off');
grid('on');

axis([0, samplRate/2.0, -25, max([max(NoiseGainL),max(AmbientGainL),max(InstrGainL), 0])]);
xlabel('Frequency, Hz','FontSize',12);
ylabel('Noise Gains, dB','FontSize',12);
title('Microphone Array Noise Gains','FontSize',14);
legend({'Total','Ambient','Instrumental'}, 'Location', 'SouthEast');


TotalOmniNoise = ComputeRMS(OmniNoise, 1, frameSize);
TotalMANoise = ComputeRMS(MANoise, 1, frameSize);
TotalNoiseGain = 20 * log10(TotalMANoise/TotalOmniNoise);
disp(sprintf('Total Noise Gain: %6.2f dB', TotalNoiseGain));

TotalAmbNoise = ComputeRMS(AmbNoise, 1, frameSize);
TotalMAAmbNoise = ComputeRMS(AmbNoise .* AmbientGain, 1, frameSize);
TotalAmbientGain = 20 * log10(TotalMAAmbNoise / TotalAmbNoise);
disp(sprintf('Total Ambient Gain: %6.2f dB', TotalAmbientGain));

TotalInstNoise = ComputeRMS(InstNoise, 1, frameSize);
TotalMAInstNoise = ComputeRMS(InstNoise .* InstrGain, 1, frameSize);
TotalInstrGain = 20 * log10(TotalMAInstNoise / TotalInstNoise);
disp(sprintf('Total Instrumental Gain: %6.2f dB', TotalInstrGain));

begFreqMeas = 300;
endFreqMeas = 3400;
begBinMeas = int32(begFreqMeas / samplRate * 2 * frameSize); 
endBinMeas = int32(endFreqMeas / samplRate * 2 * frameSize);

TotalOmniNoiseM = ComputeRMS(OmniNoise, begBinMeas, endBinMeas);
TotalMANoiseM = ComputeRMS(MANoise, begBinMeas, endBinMeas);
TotalNoiseGainM = 20 * log10(TotalMANoiseM/TotalOmniNoiseM);
disp(sprintf('Band of %d-%d Hz Noise Gain: %6.2f dB', begFreqMeas, endFreqMeas, TotalNoiseGainM));

TotalAmbNoiseM = ComputeRMS(AmbNoise, begBinMeas, endBinMeas);
TotalMAAmbNoiseM = ComputeRMS(AmbNoise .* AmbientGain, begBinMeas, endBinMeas);
TotalAmbientGainM = 20 * log10(TotalMAAmbNoiseM / TotalAmbNoiseM);
disp(sprintf('Band of %d-%d Hz Ambient Gain: %6.2f dB', begFreqMeas, endFreqMeas, TotalAmbientGainM));

TotalInstNoiseM = ComputeRMS(InstNoise, begBinMeas, endBinMeas);
TotalMAInstNoiseM = ComputeRMS(InstNoise .* InstrGain, begBinMeas, endBinMeas);
TotalInstrGainM = 20 * log10(TotalMAInstNoiseM / TotalInstNoiseM);
disp(sprintf('Band of %d-%d Hz Instrumental Gain: %6.2f dB', begFreqMeas, endFreqMeas, TotalInstrGainM));

%
%   Plot the 3D beam
%
[Fi, Theta] = meshgrid(-180:10:180);
nAngles = length(Fi);
Freq = 2500;
[val freqBin] = min(abs(Frequency - Freq)); 
WW = squeeze(Weights(curBeam, :, freqBin));

for iFi = 1:nAngles
    for iTheta = 1:nAngles
        Rho(iFi, iTheta) = ComputeArrayGain(Fi(iFi,iTheta), Theta(iFi,iTheta), Distance, Freq, MicArrDescriptor, WW);
    end
end

Fi = Fi * TORAD;
Theta = Theta * TORAD;

Rho = abs(Rho);
X = Rho .* cos(Fi) .* cos(Theta);
Y = Rho .* sin(Fi) .* cos(Theta);
Z = Rho .* sin(Theta);


fig = fig + 1;
figure(fig);
h = mesh(X,Y,Z);

light('Position',[-2,2,20])
lighting phong
material([0.4,0.6,0.5,30])
set(h,'FaceColor',[0.7 0.7 0],'BackFaceLighting','lit')
view([65,25])
axis([-1 1 -1 1 -1 1])
title('Microphone Array Directivity','FontSize',14);

end

%
%   Computes the complex array gain given azimuth, elevation, distance, frequency,
%   geometry, and weights.
%
function ArrGain = ComputeArrayGain(dir, elev, rho, freq, MicArrDescriptor, Weights)
    
    TORAD = pi / 180.0;
    soundSpeed = 342;   % m/s
    
    szMics = size(MicArrDescriptor);
    nMics = szMics(1);
    
    cX = 1;
    cY = 2;
    cZ = 3;
    aDir = 4;
    aElev = 5;
    cType = 6;
    
    sX = rho * cos(dir * TORAD) * cos(elev * TORAD);
    sY = rho * sin(dir * TORAD) * cos(elev * TORAD);
    sZ = rho * sin(elev * TORAD);
    
    OmniGain = (1/rho) * exp(-j*2.0*pi*freq * rho / soundSpeed);
    
    for iMic = 1:nMics
        dist(iMic) = sqrt((sX - MicArrDescriptor(iMic, cX))^2 + (sY - MicArrDescriptor(iMic, cY))^2 + (sZ - MicArrDescriptor(iMic, cZ))^2);
        gamma(iMic) = atan2((sY - MicArrDescriptor(iMic, cY)),(sX - MicArrDescriptor(iMic, cX)));
        gamma(iMic) = gamma(iMic) - MicArrDescriptor(iMic, aDir)*TORAD;
        distP(iMic) = sqrt((sX - MicArrDescriptor(iMic, cX))*(sX - MicArrDescriptor(iMic, cX)) + (sY - MicArrDescriptor(iMic, cY))*(sY - MicArrDescriptor(iMic, cY)));
        cappa(iMic) = atan2((sZ - MicArrDescriptor(iMic, cZ)), distP(iMic));
        cappa(iMic) = cappa(iMic) - MicArrDescriptor(iMic, aElev)*TORAD;
        CosTheta(iMic) = cos(gamma(iMic))*cos(cappa(iMic));

        Gain(iMic) = (1/dist(iMic)) * exp(-j*2.0*pi*freq * dist(iMic) / soundSpeed);
        MicDir(iMic) = MicrophoneDirectivity(dist(iMic), freq, CosTheta(iMic), MicArrDescriptor(iMic, cType));
        ASignal(iMic) = Gain(iMic) * MicDir(iMic);
    end

    OutSignal = sum(ASignal .* Weights);
    ArrGain = OutSignal / OmniGain;
end


