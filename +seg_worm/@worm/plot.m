function plot(obj)

error('Not yet finished')

% Open a big figure.
figure('OuterPosition', [50 50 1280 960]);
%fullscreen = get(0,'ScreenSize');
%figure('OuterPosition',[50 50 (fullscreen(3) - 100) (fullscreen(4) - 100)]);

% Are the head and tail flipped?
hConfidence = worm.orientation.head.confidence;
if hConfidence.head < hConfidence.tail
    worm = flipWormHead(worm);
    hConfidence = worm.orientation.head.confidence;
end

% Is the vulva on the correct side?
vConfidence = worm.orientation.vulva.confidence;
if vConfidence.vulva < vConfidence.nonVulva
    worm = flipWormVulva(worm);
    vConfidence = worm.orientation.vulva.confidence;
end

% Are we downsampling the worm?
if length(varargin) > 2
    samples = varargin{3};
else
    samples = [];
end

% When downsampling, are we interpolating the missing data or copying
% it from the original worm?
if length(varargin) > 3
    isInterp = varargin{4};
else
    isInterp = true;
end

% Downsample the worm.
if isempty(vWorm)
    if ~isempty(samples)
        origins = [0 0];
        moves = [0 0];
        pixel2MicronScale = [1 1];
        rotation = 1;
        
        % Normalize the worm.
        [vulvaContour nonVulvaContour skeleton skeletonAngles ...
            inOutTouch skeletonLength widths headArea tailArea ...
            vulvaArea nonVulvaArea] = normWorms({worm}, samples, ...
            origins, moves, pixel2MicronScale, rotation, false);
        if isInterp
            worm = [];
        end
        
        % Reconstruct the normalized worm.
        worm = norm2Worm(frame_number, vulvaContour, nonVulvaContour, ...
            skeleton, skeletonAngles, inOutTouch, skeletonLength, ...
            widths, headArea, tailArea, vulvaArea, nonVulvaArea, ...
            origins, pixel2MicronScale, rotation, worm);
    end
end

% Setup the contour, skeleton, pixels, angles, head/tail, and
% left/right sides.
contour_obj = worm.contour;
cPixels = contour_obj.pixels;
cAngles = contour_obj.angles;
wormSegSize = size(cPixels, 1) / C_WORM_SEGS;
skeleton = worm.skeleton;
sPixels = skeleton.pixels;
sAngles = worm.skeleton.angles;
cWidths = worm.skeleton.widths;
head = worm.head;
tail = worm.tail;
left = worm.left;
right = worm.right;

% Are the head and tail flipped?
if worm.orientation.head.isFlipped
    sPixels = flipud(sPixels);
    sAngles = -flipud(sAngles);
    cWidths = flipud(cWidths);
    tmp = head;
    head = tail;
    tail = tmp;
    tmp = left;
    left = right;
    right = tmp;
end

% Convert the original image to 8-bit grayscale.
if isfloat(oImg)
    oImg = uint8(round(oImg * 255));
end

% Show the original image.
hold on, subplot(2,3,1);
rgbOImg(:,:,1) = oImg;
rgbOImg(:,:,2) = oImg;
rgbOImg(:,:,3) = oImg;
imshow(rgbOImg), title('Original Image');
if ~isempty(frame_number)
    xlabel(['Frame = ' num2str(frame_number)]);
end

% Construct a binary image with the contours and skeleton overlayed.
onesImg = ones(size(img));
redImg = onesImg;
greenImg = onesImg;
blueImg = double(~img);

% Compute the unique rough and smooth contours and their intersection.
smoothCI = sub2ind(size(img), cPixels(:,1), cPixels(:,2));
roughCI = sub2ind(size(img), roughContour(:,1), roughContour(:,2));
[~, uniqueRoughI, uniqueSmoothI] = setxor(roughCI, smoothCI);
sameCI = roughCI;
sameCI(uniqueRoughI) = [];
roughCI = roughCI(uniqueRoughI);
smoothCI = smoothCI(uniqueSmoothI);

% Overlay the contour intersection.
redImg(sameCI) = 0;
greenImg(sameCI) = 0;
blueImg(sameCI) = 0;

% Overlay the rough contours.
redImg(roughCI) = 1;
greenImg(roughCI) = .7;
blueImg(roughCI) = .3;
if exist('roughIContour', 'var')
    % Ignore for now.
end

% Overlay the smooth contours and skeleton.
redImg(smoothCI) = 0;
greenImg(smoothCI) = 1;
blueImg(smoothCI) = 0;
if exist('iContour', 'var')
    % Ignore for now.
end

% Overlay the skeleton.
if ~isempty(sPixels)
    sI = sub2ind(size(img), sPixels(:,1), sPixels(:,2));
    redImg(sI) = 1;
    greenImg(sI) = 0;
    blueImg(sI) = 0;
end

% Show the binary image with the contours and skeleton overlayed.
binOImg(:,:,1) = redImg;
binOImg(:,:,2) = greenImg;
binOImg(:,:,3) = blueImg;
hold on, subplot(2,3,2);
imshow(binOImg), title(['\color{yellow}Thresholded \color{black}+ ' ...
    '(\color{orange}Rough\color{black}/\color{darkgreen}Smoothed) ' ...
    '\color{black}Contour + \color{red}Skeleton']);

% Blur the contour's local high-frequency curvature.
% Note: on a small scale, noise causes contour imperfections that shift an
% angle from its correct location. Therefore, blurring angles by averaging
% them with their neighbors can localize them better.
lfBlurSize = ceil(wormSegSize);
lfBlurWin(1:lfBlurSize) = 1 / lfBlurSize;
mcAngles = circConv(cAngles, lfBlurWin);

% Determine the min/max contour curvatures.
[mcMaxP, mcMaxI] = maxPeaksCircDist(mcAngles, lfAngleEdgeLength, ...
    cCCLengths);
[mcMinP, mcMinI] = minPeaksCircDist(mcAngles, lfAngleEdgeLength, ...
    cCCLengths);

% Determine the worm's MER (minimum enclosing rectangle).
% Note: the skeleton can exit the contour.
if ~isempty(head.bounds.contour) && ...
        ~isempty(head.bounds.skeleton) && ...
        ~isempty(tail.bounds.contour) && ...
        ~isempty(tail.bounds.skeleton) && ...
        ~isempty(left.bounds.contour) && ...
        ~isempty(left.bounds.skeleton) && ...
        ~isempty(right.bounds.contour) && ...
        ~isempty(right.bounds.skeleton)
    if isempty(sPixels)
        wMinX = min(cPixels(:,2));
        wMaxX = max(cPixels(:,2));
        wMinY = min(cPixels(:,1));
        wMaxY = max(cPixels(:,1));
    else
        wMinX = min(min(cPixels(:,2)), min(sPixels(:,2)));
        wMaxX = max(max(cPixels(:,2)), max(sPixels(:,2)));
        wMinY = min(min(cPixels(:,1)), min(sPixels(:,1)));
        wMaxY = max(max(cPixels(:,1)), max(sPixels(:,1)));
    end
    
    % Construct an image showing the head/tail, left/right sides, as well
    % as the touching/inside/outside points of the contour and skeleton.
    hRGB = [150 150 64];
    tRGB = [64 64 0];
    vRGB = [96 96 255];
    nvRGB = [0 0 224];
    cTouchRGB = [255 255 255];
    cInRGB = [255 0 0];
    cOutRGB = [0 255 0];
    sTouchRGB = [255 255 255];
    sInRGB = [0 255 0];
    sOutRGB = [255 0 0];
    sInOutRGB = [255 150 255];
    bodyImg = overlayWormTouch(oImg, worm, hRGB, 1, tRGB, 1, ...
        vRGB, 1, nvRGB, 1, cTouchRGB, 1, cInRGB, 1, cOutRGB, 1, ...
        sTouchRGB, 1, sInRGB, 1, sOutRGB, 1, sInOutRGB, 1);
    hold on, subplot(2,3,3);
    imshow(bodyImg((wMinY - 1):(wMaxY + 1),(wMinX - 1):(wMaxX + 1),:));
    if ~isempty(head.cdf) && ~isempty(tail.cdf)
        title({['Head: area=' num2str(head.area) ...
            ' cdf=[' num2str(head.cdf(1), '%.1f') ...
            num2str(head.cdf(2:end), ', %.1f') ...
            '] stdev=' num2str(head.stdev)], ...
            ['Tail: area=' num2str(tail.area) ...
            ' cdf=[' num2str(tail.cdf(1), '%.1f') ...
            num2str(tail.cdf(2:end), ', %.1f') ...
            '] stdev=' num2str(tail.stdev)]});
    end
    if ~isempty(left.cdf) && ~isempty(right.cdf)
        xlabel({['Left: area=' num2str(left.area) ...
            ' cdf=[' num2str(left.cdf(1), '%.1f') ...
            num2str(left.cdf(2:end), ', %.1f') ...
            '] stdev=' num2str(left.stdev)], ...
            ['Right: area=' num2str(right.area) ...
            ' cdf=[' num2str(right.cdf(1), '%.1f') ...
            num2str(right.cdf(2:end), ', %.1f') ...
            '] stdev=' num2str(right.stdev)]});
    end
    ylabel('Head/Tail, Left/Right, & Touch/In/Out');
    
    % Construct a pattern to identify the head.
    hImg = [1 1 1 1 1; ...
        1 1 1 1 1; ...
        1 1 1 1 1; ...
        1 1 1 1 1; ...
        1 1 1 1 1];
    [hPattern(:,1) hPattern(:,2)] = find(hImg == 1);
    hPattern(:,1) = hPattern(:,1) - ceil(size(hImg, 1) / 2);
    hPattern(:,2) = hPattern(:,2) - ceil(size(hImg, 2) / 2);
    
    % Construct a pattern to identify the vulva.
    vImg = [1 1 1 1 1; ...
        1 1 1 1 1; ...
        1 1 1 1 1; ...
        1 1 1 1 1; ...
        1 1 1 1 1];
    [vPattern(:,1) vPattern(:,2)] = find(vImg == 1);
    vPattern(:,1) = vPattern(:,1) - ceil(size(vImg, 1) / 2);
    vPattern(:,2) = vPattern(:,2) - ceil(size(vImg, 2) / 2);
    
    % Construct an image showing the contour/skeleton curvature heat map.
    blue = zeros(360, 3);
    blue(:,3) = 255;
    red = zeros(360, 3);
    red(:,1) = 255;
    cRGB = [blue(1:90,:); jet(181) * 255; red(1:90,:)]; % thermal
    sRGB = [blue(1:90,:); jet(181) * 255; red(1:90,:)]; % thermal
    sRGBNaN = [255 0 0]; % red
    hRGB = [0 255 0]; % green
    vRGB = [255 0 0]; % red
    angleImg = overlayWormAngles(oImg, worm, cRGB, sRGB, sRGBNaN, ...
        hPattern, hRGB, 1, vPattern, vRGB, 1);
    hold on, subplot(2,3,4);
    imshow(angleImg((wMinY - 1):(wMaxY + 1),(wMinX - 1):(wMaxX + 1),:));
    title({['\color{darkgreen}Head\color{black} confidence = ' ...
        num2str(hConfidence.head)], ...
        ['Tail confidence = ' num2str(hConfidence.tail)], ...
        ['Head/tail confidence = ' ...
        num2str(hConfidence.head / hConfidence.tail)]});
    xlabel({['\color{red}Vulva\color{black} confidence = ' ...
        num2str(vConfidence.vulva)], ...
        ['Non-vulva confidence = ' num2str(vConfidence.nonVulva)], ...
        ['Vulva/non-vulva confidence = ' ...
        num2str(vConfidence.vulva / vConfidence.nonVulva)]});
    ylabel('Curvature as Heat');
    
    % Show the min/max contour curvatures.
    contour_obj = worm.contour.pixels;
    hold on, text(contour_obj(mcMinI,2) - wMinX + 2, ...
        contour_obj(mcMinI,1) - wMinY + 2, '*', 'Color', 'm', ...
        'HorizontalAlignment', 'center');
    hold on, text(contour_obj(mcMinI,2) - wMinX + 2, ...
        contour_obj(mcMinI,1) - wMinY + 2, num2str(mcMinI), 'Color', 'm');
    hold on, text(contour_obj(mcMaxI,2) - wMinX + 2, ...
        contour_obj(mcMaxI,1) - wMinY + 2, '*', 'Color', 'g', ...
        'HorizontalAlignment', 'center');
    hold on, text(contour_obj(mcMaxI,2) - wMinX + 2, ...
        contour_obj(mcMaxI,1) - wMinY + 2, num2str(mcMaxI), 'Color', 'g');
end

% Biplot the contour's curvature and contour's width.
hold on, subplot(2,3,5:6);
if isempty(cWidths)
    cWidths = 0;
end
[ax h1 h2] = plotyy(1:length(cAngles), cAngles, 1:length(cWidths), ...
    cWidths);
set(h1, 'Color', 'k');
set(h2, 'Color', [.6 .6 .3]);
title(['Curvature and Width (Length = ' num2str(skeleton.length) ')']);
xlabel('Contour/Skeleton Points (Index)');
ylabel(ax(1), 'Contour/Skeleton Angle (degrees)');
ylabel(ax(2), 'Contour Width (pixels)');
%xlim(ax(1), [0 length(cAngles)]);
%xlim(ax(2), [0 length(cAngles)]);
ylim(ax(1), [-180 180]);
%maxCWidths = max(cWidths);
%ylim(ax(2), [0 maxCWidths]);
%set(ax(1), 'XTick', linspace(0, length(cAngles), 10));
set(ax(2), 'XTick', []);
set(ax(1), 'YTick', linspace(-180, 180, 13));
%set(ax(2), 'YTick', linspace(0, maxCWidths, 13));
grid on;

% Plot the contour's smoothed (min/max) curvature.
hold(ax(1), 'on'), plot(ax(1), mcAngles, 'b');
hold(ax(1), 'on'), plot(ax(1), mcMinI, mcMinP, 'm*');
hold(ax(1), 'on'), plot(ax(1), mcMaxI, mcMaxP, 'g*');

% Plot the contour's smoothed width.
if isempty(cWidths) || length(cWidths) <= 1
    mcWidths = 0;
else
    mcWidths = circConv(cWidths, lfBlurWin);
end
hold(ax(2), 'on'), plot(ax(2), mcWidths, 'c');

% Plot the skeleton's (smoothed) curvature.
% Note: we flip the skeleton angles so they're visible over the
% contour's angles.
sAngles = -sAngles;
if isempty(sAngles)
    sAngles = 0;
end
hold(ax(1), 'on'), plot(ax(1), sAngles, 'r');
msAngles = conv(sAngles, lfBlurWin, 'same');
hold(ax(1), 'on'), plot(ax(1), msAngles, 'm');

% Setup the legend.
legends = { ...
    'Contour Angles', ...
    'Avg Contour Angles', ...
    'Min Avg Contour Angles', ...
    'Max Avg Contour Angles', ...
    '- Skeleton Angles', ...
    '- Avg Skeleton Angles', ...
    'Contour Widths', ...
    'Avg Contour Widths'};
legend(legends, 'Location', 'SouthEast');

end