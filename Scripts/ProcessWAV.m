%==================================================
%
%   ProcessWAV(fileInput, fileOutput)
%
%   Processes a WAV file using overlap-add 
%
%   <fileInput>  -   input file name
%   <fileOutput> -   output file name
%
%   (c) 2007 Ivan Tashev
%
%==================================================

function ProcessWAV(fileInput, fileOutput)

%
%   Read the input data and initialization
%

[Sig, samplRate] = wavread(fileInput);
nSamples = length(Sig);
disp('Reading the input file is completed');

frameSize = 512;
frameFull = frameSize * 2;
dFrameTime = frameSize / samplRate;

prevFrame(1:frameFull) = 0.0;

nFrames = 1;

%
%   Sequentialy process the file
%
for sampleIndex = 1:frameSize:nSamples-frameFull

    dTime = nFrames * dFrameTime;
    nFrames = nFrames + 1;

    % convert to frequency domain 
    FrameSpec = ComputeFrame(Sig(sampleIndex:sampleIndex+frameFull-1));
    
    % Put frame processing code here
    OutSpec = FrameSpec;
    
    % convert the output to time doamin
    [OutSig(sampleIndex:sampleIndex+frameSize-1), prevFrame] = ReconstructFrame(OutSpec, prevFrame);

    % diagnostic output
    if (mod(nFrames, 100) == 0)
        disp(sprintf('Processed %d frames',nFrames));
    end
end

% add the last half frame
OutSig(sampleIndex:sampleIndex+frameSize-1) = prevFrame(frameSize+1:2*frameSize);
disp('Processing is completed');

%
% Write the output file on the disk
%
wavwrite(OutSig, samplRate, fileOutput);
disp('Writing the output file is completed');

