%==============================
%
%   Compute noise suppression rule
%   using various algorithms
%
%   Parameters:
%   gamma           - a posteriori SNR vector
%   xi              - a priori SNR vector
%   SuppressionType - number determining the used algorithm
%                     (see below)
%
%   (c) 2006  Ivan Tashev
%
%==============================
function Gain = SuppressionRule(gamma, xi, SuppressionType)

%
% define the suppression rule
%
SR_Nothing = 0;             % Do nothing suppression rule
SR_Wiener = 1;              % Wiener suppression rule with decision directed approach
SR_WienerApprox = 2;        % approximate Wiener suppression rule
SR_MaxLikelihood = 3;       % Maximum Likelihood suppression rule
SR_SpectralSubtraction = 4; % SpectralSubtrcation suppression rule
SR_EfraimMalah = 5;         % Efraim and Malah suppression rule
SR_EfraimMalahLog = 6;      % Efraim and Malah log suppression rule
SR_JMAP_SAE = 7;            % Joint Maximum A Posteriory Spectral Amplitude Estimator
SR_MAP_SAE = 8;             % Maximum A Posteriory Spectral Amplitude Estimator
SR_MMSE_SPE = 9;            % Minimum Mean Square Error Spectral Power Estimator

%
%   check the input data
%
if (nargin < 3)
    SuppressionType = SR_Nothing;
end
if (SuppressionType > SR_MPA_LG)
    SuppressionType = SR_Nothing;
end
    
nBins = length(xi);
if (length(gamma) ~= nBins) 
    Gain(1:nBins) = 0.0;
    return;
end

%
%   compute the suppression rule
%
nu = xi ./ (1.0 + xi) .* gamma;
eps = 1e-6;

switch SuppressionType
    case SR_Nothing                 % do nothing suppression rule
        Gain(1:nBins) = 1.0;

    case SR_Wiener                  % Wiener suppression rule
        Gain = xi ./ (1.0 + xi);

    case SR_WienerApprox            % approximate Wiener suppression rule
        Gain = max(0.0,(gamma-1.0)./(gamma+eps));

    case SR_MaxLikelihood           % Maximum Likelihood suppression rule
        Gain = 0.5 + 0.5 * sqrt(xi ./ (1.0 + xi));           

    case SR_SpectralSubtraction     % Spectral subtrcation suppression rule
        Gain = sqrt(xi ./ (1.0 + xi));
        
    case SR_EfraimMalah             % Efraim and Malah suppression rule
        Gain(1:nBins) = 1.0;        % placeholder, currently does nothing

    case SR_EfraimMalahLog          % Efraim and Malah log suppression rule
        Gain(1:nBins) = 1.0;        % placeholder, currently does nothing

    case SR_JMAP_SAE                % Joint Maximum A Posteriory Spectral Amplitude Estimator
        Gain = (xi + sqrt(xi.^2 + 2.0*(1.0+xi).*xi./(gamma+eps)))./2.0./(1.0 + xi);

    case SR_MAP_SAE                 % Maximum A Posteriory Spectral Amplitude Estimator
        Gain = (xi + sqrt(xi.^2 + (1.0+xi).*xi./(gamma+eps)))./2.0./(1.0 + xi);

    case SR_MMSE_SPE                % Minimum Mean Square Error Spectral Power Estimator
        Gain = sqrt(xi ./ (1.0 + xi) .* (1.0 + nu) ./ (gamma + eps));
end

