function angleResponse = applyAngleFilters_AllRotations(Img, SF, angle, resize)
% APPLYANGLEFILTERS_ALLROTATIONS Apply angular filter at all orientations
%
% Applies a V-shaped (angle) filter at multiple rotations to measure
% rectilinearity in an orientation-invariant way.
%
% INPUTS:
%   Img       - Input image (grayscale or RGB, double)
%   SF        - Spatial frequency (e.g., 1/32, 1/16, 1/8, 1/4)
%   angle     - Angle for the V-shaped filter (radians, 0 to pi)
%               e.g., pi/2 = 90 degree corner, pi/3 = 60 degree corner
%   resize    - Whether to resize image (1 = resize to 50%, 0 = no resize)
%
% OUTPUT:
%   angleResponse - Mean filter response at edge locations

% Filter parameters
filterSize = 300;

% Convert to grayscale if needed
if size(Img, 3) > 1
    Img = mean(Img, 3);
end

% Resize if requested
if resize == 1
    Img = imresize(Img, 0.5, 'bilinear');
end

% Edge detection
edgeImg = edge(Img, 'canny', 0.1);

% Rotation angles for orientation invariance
rotationList = 0:pi/8:(2*pi - pi/8);
nRotations = length(rotationList);

% Accumulator
angleAccum = zeros([size(Img), nRotations]);

for rot = 1:nRotations
    alpha = rotationList(rot);
    
    % Generate angle (V-shaped) filter
    [angleKernel, ~] = AngleFilter_simplified2(SF, alpha, angle, filterSize);
    angleFilt = real(angleKernel);
    
    % Apply angle filter
    filteredAngle = conv2(double(edgeImg), angleFilt, 'same');
    filteredAngle = abs(filteredAngle) / sum(abs(angleFilt(:)));
    angleAccum(:,:,rot) = abs(filteredAngle).^8;
end

% Average across orientations
angleAccum = mean(angleAccum, 3);

% Return mean response at edge locations
angleResponse = mean(angleAccum(edgeImg == 1));

end