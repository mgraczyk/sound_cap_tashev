%*******************************************
%
%   Microphone array processing program
%
%   Parameters: 
%       fileInput - string, input WAV file name
%       fileOutput - string, output WAV file name
%       direction - listeing direction, degrees, defaults 0
%
%   February 2005       Ivan Tashev
%
%*******************************************

function MicArrProcess(fileInput, fileOutput, direction)

if (nargin < 3)
    direction = 0;
end
direction = NormAngle(direction);

%=============================================
%
%   Read the input data 
%
%=============================================

disp(' ');
disp('  Microphone array processing');
disp('  (c) Ivan Tashev 2005');
disp('  ===========================');
disp(' ');
disp(sprintf('Processing: %s', fileInput));
disp(sprintf('Output File: %s', fileOutput));

% constants
TORAD = pi / 180.0;
soundSpeed = 342;   % m/s
Distance = 1.5; % m

% Read the weights
[samplRate, nMics, frameSize, nBeams, Weights, Angles] = ReadWeights('Weights.DAT');
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
nPairs = nMics - 1;

% read the wav file
[Sig, fSamplRate] = wavread(fileInput);
sSig = size(Sig);
nChans = sSig(2);
nSamples = sSig(1);
if (fSamplRate ~= samplRate)
    disp('Different sampling rate!');
    return;
end

disp('Input files reading completed');

%=============================================
%
%   Initial initialization
%
%=============================================

%-------------------------------------
% convertion to frequency domain
%-------------------------------------
frameStep = 2 * frameSize;
dFrameTime = frameSize / samplRate;
prevFrame(1:frameStep) = 0.0 + i*0.0;

%-------------------------------------
% Frequency response parameteres
%-------------------------------------
begBin     = int32( 100.0 / samplRate * 2 * frameSize);   %   zero responce below this frequency
begBinFull = int32( 200.0 / samplRate * 2 * frameSize);   %   unit gain after this frequency
endBinFull = int32(7000.0 / samplRate * 2 * frameSize);   %   unit gain up to this frequency
endBin     = int32(7500.0 / samplRate * 2 * frameSize);   %   zero response after this frequency

freqStep = samplRate / frameSize / 2;
begFreq = freqStep / 2;
Frequency(1:frameSize) = begFreq:freqStep:samplRate/2;

%-------------------------------------
%   Preprocessor
%-------------------------------------
bPreProcEnabled = 0;

%-------------------------------------
% microphone array
%-------------------------------------
bArrayEnabled = 1;
OutSpec(1:frameSize) = 0.0 + i*0.0;

%-------------------------------------
% VAD
%-------------------------------------
TrackerDataOut = SimpleVAD();
oSignal = 0;

%-------------------------------------
%   Postprocessor
%-------------------------------------
bPostProcEnabled = 0;

%-------------------------------------
% SNR estimation
%-------------------------------------
begBinMeas = int32( 300 / samplRate * 2 * frameSize); 
endBinMeas = int32(3500 / samplRate * 2 * frameSize);

nFrames = 1;
nSignalFrames = 1;
nNoiseFrames = 1;

OutNoise = 0.0;
OutSignal = 0.0;
OutSNR = 0.0;
InpNoise(1:nMics) = 0.0;
InpSignal(1:nMics) = 0.0;
InpSNR(1:nMics) = 0.0;

Ones(1:frameSize) = 1.0; 
CW = Cweighting(Frequency, Ones);
CmessW = CmessWeighting(Frequency, Ones);

%=============================================
%
%   Sequentialy process the file
%
%=============================================

%   remember the current warnings state and turn them off
s = warning('query', 'all');
warning off all

for sampleIndex = 1:frameSize:nSamples-frameStep

    dTime = nFrames * dFrameTime;
    PrevOut = OutSpec;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % convert to frequency domain 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    for iMic = 1:nMics
        FrameSpec(iMic,:) = ComputeFrame(Sig(sampleIndex:sampleIndex+frameStep-1,iMic));
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Multichannel pre-processor
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (bPreProcEnabled > 0)
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Beamformer
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (bArrayEnabled > 0)
        %apply the beamformer
        OutSpec = sum(FrameSpec .* squeeze(Weights(curBeam, :, :)));
    else
        % best microphone
        OutSpec = FrameSpec(bestMic,:);
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Clasification signal/pause
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    dLevelOut = ComputeRMS(OutSpec.*CmessW);
    TrackerDataOut = SimpleVAD(dLevelOut, dTime, TrackerDataOut);
    oSignal = TrackerDataOut.tSignal;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Post-processor
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (bPostProcEnabled > 0)
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Measure the SNR for all input channels and the output
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (nFrames > 3)
        if (oSignal > 0)
            CurOutSig = ComputeRMS(OutSpec.*CW);
            OutSignal = OutSignal + CurOutSig * CurOutSig;
            for iMic = 1:nMics
                CurInpSig = ComputeRMS(FrameSpec(iMic,:).*CW);
                InpSignal(iMic) = InpSignal(iMic) + CurInpSig * CurInpSig;
            end
            nSignalFrames = nSignalFrames + 1;
        else
            CurOutSig = ComputeRMS(OutSpec.*CW);
            OutNoise = OutNoise + CurOutSig * CurOutSig;
            for iMic = 1:nMics
                CurInpSig = ComputeRMS(FrameSpec(iMic,:).*CW);
                InpNoise(iMic) = InpNoise(iMic) + CurInpSig * CurInpSig;
            end
            nNoiseFrames = nNoiseFrames + 1;
        end
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % convert the output to time doamin
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    [OutSig(sampleIndex:sampleIndex+frameSize-1), prevFrame] = ReconstructFrame(OutSpec, prevFrame);
    
    nFrames = nFrames + 1;

    if (mod(nFrames, 500) == 0)
        disp(sprintf('Processed %d frames',nFrames));
    end
end
OutSig(sampleIndex:sampleIndex+frameSize-1) = prevFrame(frameSize+1:2*frameSize);

%   restore the warnings state
for iWarn = 1:length(s)
    warning(s(iWarn).state, s(iWarn).identifier)
end

%=============================================
%
% Display the measured SNR
%
%=============================================

disp('Generating the output signal completed');

nNoiseFrames = nNoiseFrames - 1;
nSignalFrames = nSignalFrames - 1;

for iMic = 1:nMics
    InpSignal(iMic) = sqrt(InpSignal(iMic)/nSignalFrames);
    InpNoise(iMic) = sqrt(InpNoise(iMic)/nNoiseFrames);
    InpSNR(iMic) = ComputeSNR(InpSignal(iMic), InpNoise(iMic));
    InpSignal(iMic) = 20.0 * log10(InpSignal(iMic));
    InpNoise(iMic) = 20.0 * log10(InpNoise(iMic));
    disp(sprintf(' Input %d SNR: %8.3f dBC     (%8.3f/%8.3f)', iMic, InpSNR(iMic), InpSignal(iMic), InpNoise(iMic)));
end
disp(sprintf('              %8.3f dBC Best mic', InpSNR(bestMic)));

OutSignal = sqrt(OutSignal/nSignalFrames);
OutNoise = sqrt(OutNoise/nNoiseFrames);
OutSNR = ComputeSNR(OutSignal, OutNoise);
OutSignal = 20.0*log10(OutSignal);
OutNoise = 20.0*log10(OutNoise);

disp(sprintf(' Output SNR:  %8.3f dBC     (%8.3f/%8.3f)', OutSNR, OutSignal, OutNoise));
disp(sprintf('              %8.3f dBC Improvement', OutSNR - InpSNR(bestMic)));

disp(' ');
disp(sprintf('              %8.3f dB Gain', OutSignal - InpSignal(bestMic)));
disp(' ');

%=============================================
%
% Write the output file on the disk
%
%=============================================

wavwrite(OutSig, samplRate, fileOutput);

