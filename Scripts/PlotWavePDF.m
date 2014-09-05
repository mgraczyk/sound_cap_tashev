%==================================================
%
%   PlotWavePDF(<inp_file>,<mode>)
%
%   <inp_file>  -   audio file name in WAV format
%   <mode>      -   mode to compare PDF with:
%                   G - Gaussian (default)
%                   L - Laplassian
%                   M - Gamma
%
%   (c) 2007 Ivan Tashev
%
%==================================================
function PlotWavePDF(fileName, cMode)

if (nargin < 2)
    cMode = 'G';
end
cMode = upper(cMode);

%
%   Read and normalize the WAV file
%
% read the signal
[signal1,fs] = audioread(fileName);

% remove pauses longer than 0.1 seconds
sigma1 = std(signal1);
mu1 = mean(signal1);
numSamples = length(signal1);
frame = 0.1 * fs;
sample = 1;
for sample1 = 1:frame:numSamples-frame
    if (std(signal1(sample1:sample1+frame)) > 0.05*sigma1)
        signal(sample:sample+frame) =  signal1(sample1:sample1+frame);
        sample = sample + frame;
    end
end

% normalize the signal to be zero mean and variation one
signal = signal - mean(signal);
signal = signal / std(signal);

%
%   Compute, normalize and draw the histogramm
%
points = 100;
[num, xout] = hist(signal,points);
num = num / sum(num);
plot(xout,num,'LineWidth',2);
hold on

%
%   Compute and draw the reference distribution
%
sigma = std(signal);
mu = mean(signal);

switch cMode
    case 'G'
        gauss = 1.0 / (sigma*sqrt(2*pi)) * exp(-(xout-mu).^2/(2*sigma^2));
        gauss = gauss / sum(gauss);
        plot(xout,gauss,'g','LineWidth',2);
        axis([mu-4*sigma,mu+4*sigma,0,1.2*max([max(num), max(gauss)])]);
        legend('Signal','Gauss');
       
    case 'L'
        laplass = 1.0 / (sigma * sqrt(2)) * exp(-abs(xout-mu)/(sigma*sqrt(2)));
        laplass = laplass / sum(laplass);
        plot(xout,laplass,'r','LineWidth',2);
        axis([mu-4*sigma,mu+4*sigma,0,1.2*max([max(num), max(laplass)])]);
        legend('Signal','Laplace');
        
    case 'M'
        gamma = sqrt(sqrt(3)) / (2.0*sqrt(pi*sigma)*sqrt(sqrt(2))) ./ sqrt(abs(xout)) .* exp(-sqrt(3)*abs(xout)/sqrt(2)/sigma);
        gamma = gamma / sum(gamma);
        plot(xout,gamma,'c','LineWidth',2);
        axis([mu-4*sigma,mu+4*sigma,0,1.2*max([max(num), max(gamma)])]);
        legend('Signal','Gamma');
        
    otherwise
        disp('Unsupported reference distribution');
        return;
end

%
%   Finalize the chart
%
hold off
title('Probability distribution function','FontSize',14);
xlabel('Times standard deviation','FontSize',14);
ylabel('Probability','FontSize',14);
box('on');
grid('on')
