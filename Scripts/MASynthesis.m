%*******************************************
%
%   Microphone array synthesis program
%
%   Parameters:
%
%   fileWeights    - weights file name
%                    defaults to Weights.DAT
%   fileMicArrDesc - microphone array description file name
%                    defaults to MicArrDesc.DAT
%
%   August 2005       Ivan Tashev
%
%*******************************************

function MASynthesis(fileWeights, fileMicArrDesc)

if (nargin < 1)
    fileWeights = 'Weights.DAT';
end

if (nargin < 2)
    fileMicArrDesc = 'MicArrDesc.dat';
end

%=============================================
%
%   Read the input data 
%
%=============================================

disp(' ');
disp('  Microphone array weights synthesis');
disp('  (c) Ivan Tashev 2005');
disp('  ==================================');
disp(' ');
disp(sprintf('Input from: %s', fileMicArrDesc));
disp(sprintf('Output to: %s', fileWeights));

% constants
TORAD = pi / 180.0;
soundSpeed = 342;   % m/s at 20 degrees Celsius, 1 MPa
Distance = 1.5; % m
Elev = 0.0; % elevation, degrees
eps = 1e-6;

% parameters
samplRate = 16000;
frameSize = 256;

% Read the microphone array geometry
MicArrDescriptor = load(fileMicArrDesc);
cX = 1;
cY = 2;
cZ = 3;
aDir = 4;
aElev = 5;
cType = 6;
szMics = size(MicArrDescriptor);
nMics = szMics(1);

% Frequency response parameteres
begBin = int32(100 / samplRate * 2 * frameSize);        %   zero responce below this frequency
begBinFull = int32(200 / samplRate * 2 * frameSize);    %   unit gain after this frequency
endBinFull = int32(7000 / samplRate * 2 * frameSize);   %   unit gain up to this frequency
endBin = int32(7500 / samplRate * 2 * frameSize);       %   zero response after this frequency

FreqResponse(1:begBin) = 0;
for iBin = begBin:begBinFull
    FreqResponse(iBin) = 0.5 * (1.0 + cos(double(iBin - begBinFull)/double(begBinFull - begBin)*pi));
end
FreqResponse(begBinFull:endBinFull) = 1;
for iBin = endBinFull:endBin
    FreqResponse(iBin) = 0.5 * (1.0 + cos(double(iBin - endBinFull)/double(endBin - endBinFull)*pi));
end
FreqResponse(endBin:frameSize) = 0;

freqStep = samplRate / frameSize / 2;
begFreq = 0.0;          % for FFT
Frequency(1:frameSize) = begFreq:freqStep:samplRate/2-freqStep/2;

% determine the work range
stepAngle = 10.0;    % degreees
begAngle = -50;
endAngle = 50;

Angles = begAngle:stepAngle:endAngle;
nBeams = length(Angles);

Weights(1:nBeams,1:nMics,1:frameSize) = 0.0 + i * 0.0;

disp('Initializing completed');

%=============================================
%
%   In-line function for DS weights computation
%
%=============================================
function Weight = DSBeamDesign(iBin)
    Weight(1:nMics) = 0.0;
    
    for iMic = 1:nMics
        Weight(iMic) = OmniGain / InpSignal(iMic);
        Weight(iMic) = Weight(iMic) / abs(Weight(iMic));
    end

    OutSignal = sum(Weight .* InpSignal);
    NormGain = OmniGain / OutSignal;
    Weight = NormGain * Weight;
end

%=============================================
%
%   In-line function for BM weights computation
%
%=============================================
function Weight = BMBeamDesign(iBin)
    Weight(1:nMics) = 0.0;
    Weight(bestMic) = OmniGain / InpSignal(bestMic);
    
    OutSignal = sum(Weight .* InpSignal);
    NormGain = OmniGain / OutSignal;
    Weight = NormGain * Weight;
end

%=============================================
%
%   In-line function for initial point computation
%
%=============================================
function Weight = PreBeamDesign(iBin)
    Weight = DSBeamDesign(iBin);
    WeightBM = BMBeamDesign(iBin);

    NoiseGainDS = ComputeNoiseGain(Weight, MicArrDescriptor, Freq);
    NoiseGainBM = ComputeNoiseGain(WeightBM, MicArrDescriptor, Freq);
    if (NoiseGainDS > NoiseGainBM)
        Weight = WeightBM;
        disp(sprintf('     Replaced with best mic at %7.1f Hz',Freq));
    end
end

%=============================================
%
%   In-line function for one bin design
%
%=============================================
function Weight = BeamDesign(iBin, Weight0)
    % Initialization
    Freq = Frequency(iBin);
    OmniGain = (1.0/Distance) * exp(-j*2.0*pi*Freq * Distance / soundSpeed);

    for iMic = 1:nMics
        Gain(iMic) = (1/dist(iMic)) * exp(-j*2.0*pi*Freq * dist(iMic) / soundSpeed);
        MicDir(iMic) = MicrophoneDirectivity(dist(iMic), Freq, CosTheta(iMic), MicArrDescriptor(iMic, cType));
        InpSignal(iMic) = Gain(iMic) * MicDir(iMic);
    end

    Weight = PreBeamDesign(iBin);
    NoiseGain = ComputeNoiseGain(Weight, MicArrDescriptor, Freq);
    
    disp(sprintf('  Freq = %7.1f Hz, NG = %6.4f', Freq, NoiseGain));
end

%=============================================
%
%   Estimation loop
%
%=============================================

%   remember the current warnings state and turn them off
s = warning('query', 'all');
warning off all

for iBeam = 1:nBeams

    %
    %   Beam related computations
    %
    dir = Angles(iBeam);
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
    
    bestMic = BestMicrophone(MicArrDescriptor, dir, Distance, Elev);
    
    disp(sprintf('Design beam at %d degrees [%d]', dir, bestMic));
    
    firstBin = 1;
    
    Freq = 0.0;
    OmniGain = 0.0;
    InpSignal(1:nMics) = 0.0;

    %
    %   Per frequency bin design
    %
    for iBin = begBin:endBin
        Weights(iBeam, :, iBin) = BeamDesign(iBin, squeeze(Weights(iBeam, :, iBin-1)));
    end
    
    % normalize
    for iBin = begBin:endBin
        Weights(iBeam, :, iBin) = FreqResponse(iBin) .* squeeze(Weights(iBeam, :, iBin));
    end
   
end

%   restore the warnings state
for iWarn = 1:length(s)
    warning(s(iWarn).state, s(iWarn).identifier)
end

disp('Computing the responce completed');

WriteWeights(fileWeights, samplRate, nMics, frameSize, Weights, Angles);

end
