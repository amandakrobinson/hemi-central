function [SpaceKernel, FreKernel] = AngleFilter_simplified(SF, alpha, theta, filterSize)
% ANGLEFILTER_SIMPLIFIED Generate banana wavelet kernels for rectilinearity measurement
%
% At alpha=0, theta=0: horizontal Gabor bar with horizontal fringes.
% Increasing theta opens the V toward a corner.
% alpha rotates the whole filter in the image.
%
% INPUTS:
%   SF         - Spatial frequency (e.g., 1/32, 1/16, 1/8, 1/4)
%   alpha      - Rotation angle of the filter in the image (radians)
%   theta      - Corner angle: 0 = straight line, pi/2 (90°) = right angle corner,
%                pi*5/6 (150°) = obtuse corner, approaching pi = very acute corner
%   filterSize - Half-size of filter kernel (e.g., 300 creates 600x600 kernel)
%
% OUTPUTS:
%   SpaceKernel - Complex-valued spatial domain filter
%   FreKernel   - Frequency domain representation
%
% Based on banana wavelet method from:
%   Kruger, N. - banana wavelets
%   Nasr et al., 2014 - rectilinearity in scene perception

if nargin < 4
    filterSize = ceil(8 / SF);
end

% Wave vector magnitude
k = 2 * pi * SF;

% Bending parameter
% theta=0 -> straight (bA=0), theta=pi/2 -> right angle, theta->pi -> very acute
theta_clamped  = max(min(theta, pi - 0.01), 0.01);
theta_internal = pi - theta_clamped;
bA = 1 / tan((pi/2) - (theta_internal/2));

% Coordinate grid
[xGrid, yGrid] = meshgrid(-filterSize:filterSize-1, -filterSize:filterSize-1);

% Rotate by alpha
xRot = cos(alpha)*xGrid + sin(alpha)*yGrid;
yRot = cos(alpha)*yGrid - sin(alpha)*xGrid;

% --- Two-arm Gaussian envelope ---
% Each arm lies along direction (±1, bA) in (xRot, yRot) space.
% armAngle: angle each arm makes with the xRot axis
armAngle = atan2(bA, 1);

% Scale with filterSize for consistent size across all SF
wavelength  = 1 / SF;
sigmaAlong  = wavelength * 2.0;   % length along arm
sigmaAcross = wavelength * 0.5;   % thickness across arm

% Arm 1: xRot <= 0, pointing in direction (-cos(armAngle), sin(armAngle))
along1  = -xRot * cos(armAngle) + yRot * sin(armAngle);
across1 =  xRot * sin(armAngle) + yRot * cos(armAngle);
gauss1  = exp(-0.5 * ((along1/sigmaAlong).^2 + (across1/sigmaAcross).^2)) .* (xRot <= 0);

% Arm 2: xRot > 0, pointing in direction (cos(armAngle), sin(armAngle))
along2  =  xRot * cos(armAngle) + yRot * sin(armAngle);
across2 = -xRot * sin(armAngle) + yRot * cos(armAngle);
gauss2  = exp(-0.5 * ((along2/sigmaAlong).^2 + (across2/sigmaAcross).^2)) .* (xRot > 0);

% Windowing function
windowSigma = filterSize * 0.4;
windowFunc  = exp(-(xGrid.^2 + yGrid.^2) / (2 * windowSigma^2));

% Complex carrier — uses across coordinate so fringes run parallel to each arm
% At alpha=0, theta=0: across = yRot, so fringes are horizontal
realPart = (gauss1 .* cos(k * across1) + gauss2 .* cos(k * across2)) .* windowFunc;
imagPart = (gauss1 .* sin(k * across1) + gauss2 .* sin(k * across2)) .* windowFunc;

% DC-free
realPart = realPart - mean(realPart(:));
imagPart = imagPart - mean(imagPart(:));

% Unit energy normalisation
energy = sqrt(sum(realPart(:).^2) + sum(imagPart(:).^2));
if energy > 0
    realPart = realPart / energy;
    imagPart = imagPart / energy;
end

SpaceKernel = complex(realPart, imagPart);
FreKernel   = fft2(SpaceKernel);

end