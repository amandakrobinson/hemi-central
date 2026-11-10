function make_rdm_models()

% for each stimulus, calculate:
% -- mean luminance
% -- spatial frequency values
% -- rectilinearity
% -- colour
% then construct RDMs for each feature using euclidean distance

% also make category models for concept and object/word

load('results/imageorder.mat','imnames')

[edges,cielab,lum,sflow,sfhigh] = deal([]);

for i = 1:length(imnames)
    
    % counter
    imname = dir(sprintf('stimuli_rendered/stim%03i*',i));
    fprintf('Calculating values for %s\n',imname.name);
    
    % read image
    ii = imread(sprintf('stimuli_rendered/%s',imname.name));
    if size(ii,3) == 1
        ii = cat(3,ii,ii,ii);
    end
    iig = rgb2gray(im2double(ii));
    
    % outline
    edges(:,:,i)=edge(iig,'canny',.1);
    
    % luminance
    lum(i) = mean2(iig);
    
    % colour
    ii = rgb2lab(ii); % convert to CIElab

    % bin a* and b* channels into a 2D histogram (channel 1 is luminance -
    % don't need that here since we have a separate luminance model)
    a_chan = ii(:,:,2);
    b_chan = ii(:,:,3);
    edges_chans = linspace(-128, 128, 33); % 32 bins for each channel
    hist_ab = histcounts2(a_chan(:), b_chan(:), edges_chans, edges_chans);
    hist_ab = hist_ab(:) / sum(hist_ab(:)); % normalise
    cielab(i,:) = hist_ab;

    %% spatial frequency

    % calculate radial SF - direction invariant
    % images were 128 x 128 and presented at 3.24 dova
    % display pixels per degree: 128 / 3.24 = 39.5 pixels/degree
    % display Nyquist frequency: 39.5 / 2 = 19.75 cycles/degree (64 cycles per image)
    % frequency per index: 19.75/64 = 0.309 cycles/degree

    % radial bins: index k contains radii (k-2, k-1] cycles/image; index 1 = DC
    % Low SF  (<3 cpd):        indices 2:10  = 0.31-2.78 cpd
    % High SF (6 cpd-Nyquist): indices 22:65 = 6.17-19.75 cpd


    % Mean-subtract image before FFT to remove DC component
    % This ensures low SF captures contrast energy at low frequencies
    % rather than overall brightness (which is captured by the lum model)
    iig_ms = iig - mean2(iig);


    % Perform 2D FFTs
    fftOriginal = fft2(iig_ms);
    
    % Move center from (1,1) to the middle of the matrix.
    shiftedFFT = fftshift(fftOriginal);
    shiftedFFTMagnitude = abs(shiftedFFT);

    % Get the average radial profile
    [rows, columns] = size(iig);
    midRow = rows/2+1;
    midCol = columns/2+1;
    maxRadius = ceil(sqrt(midRow^2 + midCol^2));
    radialProfile = zeros(maxRadius, 1);
    count = zeros(maxRadius, 1);
    for col = 1 : columns
        for row = 1 : rows
            radius = sqrt((row - midRow) ^ 2 + (col - midCol) ^ 2); % how far is this pixel from the centre
            thisIndex = ceil(radius) + 1;
            radialProfile(thisIndex) = radialProfile(thisIndex) + shiftedFFTMagnitude(row, col);
            count(thisIndex) = count(thisIndex) + 1;
        end
    end

    % Get average
    radialProfile = radialProfile ./ count;
    
    % Low SF (0-3 cpd): start from index 2 to skip DC component
    % Index range: 2 to floor(3 / 0.309)+1 = 2 to 10
    sflow(i,:) = radialProfile(2:10);

    % High SF (6+ cpd):
    % Index range: ceil(6 / 0.309)+2 to ceil(19.75/0.309)+1 = 22 to 65
    sfhigh(i,:) = radialProfile(22:65);

end

%% now get rectilinearity
load('results/imagerect.mat','normData')

rect = reshape(normData,36,[]); % already in that form

%% make image model RDMs

% low-level models
m_cielab = squareform(pdist(cielab,'euclidean'));
m_lum = squareform(pdist(lum','euclidean'));
m_sflow = squareform(pdist(sflow,'euclidean'));
m_sfhigh = squareform(pdist(sfhigh,'euclidean'));
m_rect = squareform(pdist(rect,'euclidean'));
m_edge = squareform(pdist(reshape(edges,[],36)','jaccard'));

%% make category model RDMs

% higher level category model
m_cat = ones(36,36);
for c = 1:6
    init = (c-1)*4+1;
    widx = 25+(c-1)*2;
    idx = [init:(init+3) widx widx+1];
    m_cat(idx,idx) = 0;
end

% image vs word model
m_imw = ones(36,36);
m_imw(1:24,1:24) = 0;
m_imw(25:36,25:36) = 0;

save('results/RSAmodels.mat','m_cielab','m_sflow','m_sfhigh','m_rect','m_lum', ...
    'm_cat','m_imw','m_edge');

   