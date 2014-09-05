%*******************************************
%
%   Microphone array processing program
%
%   Parameters: 
%       fileInput - string, input WAV file name
%       fileOutput - string, outnput WAV file name
%       direction - listeing directin, degrees, defaults 0
%
%   February 2005       Ivan Tashev
%
%*******************************************

function MicArrSSL(fileInput)

%=============================================
%
%   Read the input data 
%
%=============================================

disp(' ');
disp('  Microphone array sound source localization');
disp('  (c) Ivan Tashev 2008');
disp('  ==========================================');
disp(' ');
disp(sprintf('Processing: %s', fileInput));

% constants
TORAD = pi / 180.0;
soundSpeed = 342;   % m/s
Distance = 1.5; % m

% Read the microphone array geometry
MicArrDescriptor = load('MicArrDesc.dat');
cX = 1;
cY = 2;
cZ = 3;
aDir = 4;
aElev = 5;
cType = 6;
szMics = size(MicArrDescriptor);
nMics = szMics(1);
nPairs = nMics*(nMics - 1)/2;

% read the wav file
[Sig, samplRate] = wavread(fileInput);
sSig = size(Sig);
nChans = sSig(2);
nSamples = sSig(1);
if (nChans < nMics)
    disp('Insufficient number of channels!');
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
frameSize = 256;
frameStep = 2 * frameSize;
dFrameTime = frameSize / samplRate;
prevFrame(1:frameStep) = 0.0 + i*0.0;
nFrames = 0;

%-------------------------------------
% Frequency response parameteres
%-------------------------------------
begBin     = int32( 200.0 / samplRate * 2 * frameSize); 
endBin     = int32(7000.0 / samplRate * 2 * frameSize); 

freqStep = samplRate / frameSize / 2;
begFreq = freqStep / 2;
Frequency(1:frameSize) = begFreq:freqStep:samplRate/2;

%-------------------------------------
% VAD
%-------------------------------------
TrackerDataOut = SimpleVAD();
oSignal = 0;

%-------------------------------------
% SNR estimation
%-------------------------------------
Ones(1:frameSize) = 1.0; 
CW = Cweighting(Frequency, Ones);
CmessW = CmessWeighting(Frequency, Ones);

disp('Initialization completed');

%-------------------------------------
% Sound source localization
%-------------------------------------
SSLframe = 0;
NoiseVar(1:nMics,1:frameSize) = 0.0;

hypAngles = -90:5:90;
numHypAngles = length(hypAngles);

nPair = 0;
for iMic1 = 1:nMics-1
    for iMic2 = iMic1+1:nMics
        nPair = nPair + 1;
        dist(nPair) = sqrt(sum((MicArrDescriptor(iMic1,cX:cY)-MicArrDescriptor(iMic2,cX:cY)).^2));
    end
end
d = max(dist);
maxTaps = fix(d/soundSpeed*samplRate);

% single pair
TDE_1Pair = 0;

% TDE algorithms
TDE_Interp = 1;
TDE_Estim = 2;

% SRP algorithms
SRP_BF = 3;
SRP_PHAT = 4;
SRP_MLR = 5;
SRP_ML = 6;
SRP_MUSIC = 7;
SRP_PHAT_Beam = 10;

SSLalgorithm = TDE_1Pair;

K = 10;
XX(1:K, 1:nMics,1:frameSize) = 0.0;
PP(1:numHypAngles,1:frameSize) = 0.0;

disp('SSL Initialization completed');

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
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % convert to frequency domain 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    for iMic = 1:nMics
        FrameSpec(iMic,:) = ComputeFrame(Sig(sampleIndex:sampleIndex+frameStep-1,iMic));
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Clasification signal/pause
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    dLevelOut = mean((mean(abs(FrameSpec)).*CmessW).^2);
    TrackerDataOut = SimpleVAD(dLevelOut, dTime, TrackerDataOut);
    oSignal = TrackerDataOut.tSignal;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Sound source localizer
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %keep the last K frames
    if ((SSLalgorithm == SRP_BF) || (SSLalgorithm == SRP_PHAT) || (SSLalgorithm == SRP_MLR) ||(SSLalgorithm == SRP_ML) ||(SSLalgorithm == SRP_MUSIC))
        XX = circshift(XX,1);
        XX(1,:,:) = FrameSpec;
    end
    
    if (oSignal > 0)
        SSLframe = SSLframe + 1;
        SSLtime(SSLframe) = dTime;
        
        % compute the spectral matrix
        for iBin = 1:frameSize
            X = squeeze(XX(:,:,iBin)).';
            S = X*X';
            SS(iBin,:,:) = S;
        end
        
        %
        % compute the spatial response function
        %
        P(1:numHypAngles) = 0.0;
        switch (SSLalgorithm)
            
            %
            %  Time delay estimation with GCC, use channels 1 and 2
            % 
            case TDE_1Pair
                P = ComputeGCC(squeeze(FrameSpec(1,:)), squeeze(FrameSpec(2,:)), squeeze(NoiseVar(1,:)), squeeze(NoiseVar(2,:)), begBin, endBin, frameSize, maxTaps);

            %
            %  Time delay estimation based algorithms
            % 
            % Combining the pairs by summing the GCCs after linear interpolation
            case TDE_Interp

            % Computing the GCC in one desired point - equivalent to DS beamsteering
            case TDE_Estim
                
            %
            %   Steered response power methods
            %
            % traditional steered response power algorithm
            case SRP_BF

            % PHAT weighted steered response power algorithm
            case SRP_PHAT
                
            % MLR weighted steered response power algorithm
            case SRP_MLR

            % PHAT weighted MVDR steered response power algorithm (actually MPDR)
            case SRP_ML

            % MUSIC algorithm
            case SRP_MUSIC
                
        end

        PPP(SSLframe,:) = P;
    else
        % compute the noise variation
        NoiseVar = 0.98 * NoiseVar + 0.02 * abs(FrameSpec).^2;
    end
    
    nFrames = nFrames + 1;
    if (mod(nFrames, 500) == 0)
        disp(sprintf('Processed %d frames',nFrames));
    end
end

%   restore the warnings state
for iWarn = 1:length(s)
    warning(s(iWarn).state, s(iWarn).identifier)
end

%=============================================
%
% Display the measured SNR
%
%=============================================
disp('Processing the input file completed');

%=============================================
%
% Write the output file on the disk
%
%=============================================

save('PPP.dat','PPP','-ascii');
save('SSLTime.dat','SSLtime','-ascii');

end

%
%   Compute Generalized Crocc-Correlation function with various weightings
%
function GCCpart = ComputeGCC(Sig1, Sig2, Noise1, Noise2, begBin, endBin, frameSize, maxTaps)

    % Set TDE weigting type
    GCC_NO = 1;
    GCC_Roth = 2;
    GCC_SCOT = 3;
    GCC_PHAT = 4;
    GCC_Eckhart = 5;
    GCC_ML = 6; 
    GCC_MLA = 7;
    GCC_MLR = 8;

    GCCweigting = GCC_PHAT;

    % compute the auto and cross correlation function spectra
    Gx1x2 = Sig1.*conj(Sig2);
    Gx1x1 = Sig1.*conj(Sig1);
    Gx2x2 = Sig2.*conj(Sig2);

    % weight the CC
    switch (GCCweigting)
        case GCC_NO
            Psi(1:frameSize) = 1.0;

        case GCC_Roth
            Psi = 1.0 ./ (Gx1x1 + eps); 

        case GCC_SCOT
            Psi = 1.0 ./ (sqrt(Gx1x1 .* Gx2x2)+ eps);

        case GCC_PHAT
            Psi = 1.0 ./ (abs(Gx1x2) + eps);

        case GCC_Eckhart
            Psi = abs(Gx1x2) ./ (max(eps,(Gx1x1 - abs(Gx1x2))) .* max(eps,(Gx2x2 - abs(Gx1x2)))); 

        case GCC_ML
            G12 = Gx1x2 ./ (sqrt(Gx1x1 .* Gx2x2)+ eps);
            absG12sq = abs(G12).^2;
            Psi = absG12sq ./ ((abs(Gx1x2).*(1-absG12sq))+eps);

        case GCC_MLA
            Psi = Sig1.*Sig2 ./ ((Noise1.*(Sig2.^2) + Noise2.*(Sig1.^2))+eps);

        case GCC_MLR
            q = 0.25;
            Psi = 1.0 ./ (q * abs(Gx1x2) + (1-q)*0.5*(Noise1+Noise2));
    end

    % compute the CC in time domain
    Psi(1:begBin) = 0.0;
    Psi(endBin:frameSize) = 0.0;
    GCC = Psi .* Gx1x2;
    GCC(frameSize+1) = 0.0 + i * 0.0;
    GCC(frameSize+2:2*frameSize) = conj(GCC(frameSize:-1:2));
    RCC = ifftshift(ifft(GCC));

    % select the center part
    begIndex = frameSize-maxTaps;
    endIndex = frameSize+maxTaps+2;

    GCCpart = RCC(begIndex:endIndex);

end

