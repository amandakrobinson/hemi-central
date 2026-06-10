
%% Calculate rectilinearity score for images
% Convolve images with angular filters at different spatial frequencies 
% and orientations to quantify rectilinearity.
%
% Based on: Nasr, S., Echavarria, C.E. and Tootell, R.B., 2014. 
% Thinking outside the box: rectilinear shapes selectively activate 
% scene-selective cortex. Journal of Neuroscience, 34(20), pp.6721-6735.


%% Image details
picfolder = 'stimuli_rendered';
imfiles = dir(sprintf('%s/*.png', picfolder));
nPics = length(imfiles);
ims = {imfiles.name};

%% Filter parameters
% Spatial frequencies
SFlist = [1/24  1/12 1/6];

% Angles for V-shaped filters (radians)
angleList = [0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6];  % ~0 to 150 in 30 deg steps
angleLabels = {'0°', '30°', '60°', '90°', '120°', '150°'};

nSF = length(SFlist);
nAngles = length(angleList);

%% Visualise the filters
figure(3); clf;
set(gcf, 'Position', [1625 -43 1800 1000]);

filterSize = 300;
alphaA = 0;%3*pi/2;

for s = 1:nSF
    SF = SFlist(s);
    for a = 1:nAngles
        angle = angleList(a);
        [angleK, ~] = AngleFilter_simplified2(SF, alphaA, angle, filterSize);
        subplot(nSF, nAngles, (s-1)*nAngles + a);
        imagesc(real(angleK));
        axis image; axis off;
        title(sprintf('Angle: %.0f°', rad2deg(angle)));
    end
end
sgtitle('Angular filters', 'FontSize', 13);
saveas(gcf, 'figures/filters.png');

%% Get rectilinearity scores per image
outputMatrix = zeros(nPics, nSF, nAngles);

for pic = 1:nPics
    ImgName = [picfolder '/' ims{pic}];
    fprintf('Analysing rectilinearity in image %s\n', ims{pic});
    im = imread(ImgName);
    if size(im, 3) == 1
        im = repmat(im, 1, 1, 3);
    end
    Img = rgb2gray(im2double(im));
    
    for s = 1:nSF
        SF = SFlist(s);
        for a = 1:nAngles
            angle = angleList(a);
            outputMatrix(pic, s, a) = applyAngleFilters_AllRotations(Img, SF, angle, 0);
        end
    end
end

%% Normalize rectilinearity scores
normData = zeros(size(outputMatrix));

for a = 1:nAngles
    for s = 1:nSF
        temp = outputMatrix(:, s, a);
        minVal = min(temp(:));
        maxVal = max(temp(:));
        if maxVal > minVal
            normData(:, s, a) = (temp - minVal) / (maxVal - minVal);
        else
            normData(:, s, a) = 0;
        end
    end
end

finalData = squeeze(mean(mean(normData, 2), 3));

%% Plot rectilinearity details
figure(1); clf;

subplot(2, 2, 1);
imagesc(finalData);
title('Overall rect index');

subplot(2, 2, 2);
imagesc(mean(normData, 3));
title('Norm rect by spatial scale');

for a = 1:nAngles
    subplot(2, nAngles, nAngles + a);
    imagesc(normData(:, :, a));
    title(sprintf('Angle %s', angleLabels{a}));
    ylabel('SF');
end

saveas(gcf, 'figures/image_rectilinearity_values.png');

%% Plot images by rank
[val, ind] = sort(finalData, 'ascend');

f = figure(2); clf;
f.Position = [0 0 1917 863];

for p = 1:nPics
    subplot(4, 9, p);
    Img = imread([picfolder '/', ims{ind(p)}]);
    imshow(Img);
    title(sprintf('%.02g', val(p)));
end

saveas(gcf, 'figures/image_rectilinearity.png');

%% Save data
save('imagerect.mat', 'normData', 'outputMatrix', 'finalData', ...
     'SFlist', 'angleList', 'angleLabels', 'ims');