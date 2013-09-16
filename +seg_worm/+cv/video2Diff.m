function video2Diff(video_file_path, varargin)
%video2Diff Differentiate a video frame by frame. 
%
%   The difference is defined as the variance of pixel differences between
%   subsequent frames. Note 1: video frame indexing begins at zero;
%   therefore, video differentiation indexing begins at 1. Note 2: dropped
%   frames are labelled as NaN. Frames separated by dropped frames are
%   differentiated as subsequent. Therefore, these frames can be recognized
%   as any value following a NaN.
%
%
%   seg_worm.cv.video2Diff(videoFile, *diffVideoFile, *progressFunc, *progressState, *remMode, *dilatePixels)
%
%   Inputs:
%       videoFile     - the name of the video to differentiate

%       diffVideoFile - the name for the video of the frame differences;
%                       if empty, we don't output a video
%       progressFunc  - the function to call as an update on our progress.
%                       Video differentiation is slow. The progress
%                       function, if not empty, is called every frame. This
%                       function can keep a GUI updated on the progress of
%                       video differentiation. The progress function must
%                       have the following signature:
%
%                       progressFunc(state, frames, frame, image,
%                                    diffImage, variance)
%
%                       state     = persistent state data for the function
%                       frames    = the total number of frames
%                       frame     = the frame number
%                       image     = the video image
%                       diffImage = the frame-difference image
%                       variance  = the frame-difference variance
%       progressState - the persistent state data for the progressFunc
%       remMode       - the worm removal mode. Occasionally fast worm
%                       movements can be confused with stage movements.
%                       Removing the worm eliminates its contribution to
%                       frame differences (unfortunately, this computation
%                       is temporally expensive). The modes are:
%
%                       0 = don't remoe the worm; the default
%                       1 = remove everything darker than the Otsu threshold
%                       2 = remove the largest Otsu-thresholded connected component
%
%       dilatePixels  - the number of neighboring pixels to use for
%                       dilating the lagest connected component (used in
%                       conjunction with worm removal); if empty, the
%                       largest connected component isn't dilated
%
%   Old Inputs
%   -----------------------------
%   diffFile      - the name for the file containing the frames/second
%                       and the video differences per subsequent frames
%                       (the difference is defined as the variance of
%                       pixel differences between subsequent frames).
%   JAH NOTE: I removed diffFile and made this an automatic name
%   I think it was optional because you can have all of the different
%   options. I might change this back later ...
%       
%
% See also:
%   FINDSTAGEMOVEMENT
%   SIMPLEPROGRESSUPDATE
%
%   JAH NOTES: This function is important for finding stageMovements
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.


% Are we creating a video of the frame differences?
diffVideoFile = [];
if ~isempty(varargin)
    diffVideoFile = varargin{1};
end

% Is there a progress function?
progressFunc = [];
if length(varargin) > 1
    progressFunc = varargin{2};
end

% Is there a state for the progress function?
progressState = [];
if length(varargin) > 2
    progressState = varargin{3};
end

% Determine the thresholded-object removal mode.
remMode = 0;
if length(varargin) > 3
    remMode = varargin{4};
end

% Create the dilation disk.
dilateDisk = [];
if length(varargin) > 4 && varargin{5} > 0
    dilateDisk = strel('disk', varargin{5});
end

%--------------------------------------------------------------------------

% Open the video and get the seconds / frame.
vr = videoReader(video_file_path,true);

fps    = vr.fps;
n_frames = vr.n_frames;

% Check the frame rate.
%--------------------------------------------------
minFPS = .1;
maxFPS = 100;
if fps < minFPS || fps > maxFPS
    warning('video2Diff:WeirdFPS', [videoFile ' was recorded at ' ...
        num2str(fps) ' frames/second. An unusual frame rate']);
end

file_manager = seg_worm.file_manager(video_file_path);

% Get the vignette.
%--------------------------------------------------
if remMode > 0 
   vignette     = seg_worm.vignette.create(file_manager,vr);
else
   vignette = [];
end

%{
%Are we making a video of the differences?
%-------------------------------------------------------
if ~isempty(diffVideoFile)
    
    % Construct the video file.
    diffVideoFile = strrep(videoFile, '.avi', '_diff.avi');
    vw = videoWriter(diffVideoFile, 'fps', fps, 'plugin', 'DirectShow');
    
    % Initialize the video information.
    colorOffset = 16;
    noImg = uint8(zeros(size(img)));
    
    % Construct an empty RGB image.
    noRGBImg(:,:,1) = noImg;
    noRGBImg(:,:,2) = noImg;
    noRGBImg(:,:,3) = noImg;
    prevDiffRGBImg = noRGBImg;
end
%}

img = vr.getFrame(1);
img = single(helper__handleRemoveMode(img,remMode,vignette,dilateDisk));
last_valid_image = img;

frameDiffs = NaN(1,n_frames-1,'single');

for cur_frame = 2:n_frames
    
    img = vr.getFrame(cur_frame);
    %NOTE: img might be empty ...
    
    img = single(helper__handleRemoveMode(img,remMode,vignette,dilateDisk));
    
    if ~isempty(img)
        diffImg = img - last_valid_image;
        
        %NOTE: We shouldn't get NaN values if we are not removing anything
        %...
        if remMode > 0
            mask    = ~isnan(diffImg(:)); %nanvar() slow, we'll remove nan ourselves
            x = diffImg(mask);
        else
            x = diffImg(:); 
        end
        
        %From var() - we remove var() since it is a suprisingly slow call
        %-----------------------------------------------------------
        %We assume:
        %1) We always have at least one valid value
        %2) We are operating along a single dimension
        n     = single(length(x));
        denom = n - 1;
        
        %   var() code
        %   ---------------------------------------------------------------
        %   xbar = sum(x, dim) ./ n;
        %   x = bsxfun(@minus, x, xbar);
        %   y = sum(abs(x).^2, dim) ./ denom; % abs guarantees a real result
        %------------------------------------------------------------------
        %   x = bsxfun(@minus,x,sum(x)./n) <- combining top 2 lines
        %   -> don't need bsxfun - 1d
        %   -> don't need abs - real values
        
        frameDiffs(cur_frame-1) = sum((x - (sum(x)./n)).^2)./denom;

        if remMode > 0
            %???? - why do this? - just because ??????
            frameDiffs(cur_frame-1) = frameDiffs(cur_frame-1) .^ 2;
        end
        last_valid_image = img; 
    end
    
    % Update the progress.
    if ~isempty(progressFunc)
        index = cur_frame - 1;
        progressFunc(progressState, n_frames, index, img, diffImg, frameDiffs(index)); %#ok<NOEFF>
    end
    
    %{
    % Are we making a video of the differences?
    if ~isempty(diffVideoFile)
        
        % Label the negative changes in red.
        rImg = noImg;
        rDiffs = diffImg < 0;
        rImg(rDiffs) = uint8(abs(diffImg(rDiffs)) + colorOffset);
        
        % Label the positive changes in green.
        gImg         = noImg;
        gDiffs       = diffImg > 0;
        gImg(gDiffs) = uint8(abs(diffImg(gDiffs)) + colorOffset);
        
        % Construct the RGB image.
        diffRGBImg(:,:,1) = rImg;
        diffRGBImg(:,:,2) = gImg;
        diffRGBImg(:,:,3) = noImg;
        
        % Add the difference frame(s).
        for k = 1:droppedFrames
            addframe(vw, prevDiffRGBImg);
        end
        addframe(vw, diffRGBImg);
        prevDiffRGBImg = diffRGBImg;
    end
    %}
end

% Clean up.
close(vr);

if ~isempty(diffVideoFile)
    close(vw);
end

% Save the video differentiation.
save(file_manager.diff_file, 'fps', 'frameDiffs');


end

function img = helper__handleRemoveMode(img,remMode,vignette,dilateDisk)

if isempty(img)
   return 
end


% Remove dark objects from the image.
if remMode > 0
    
    % Correct the vignette.
    if ~isempty(vignette)
        img = vignette.apply(img);
    end
    
    % Binarize the image.
    bImg = otsuImg(img, true);
    bImg = ~bImg;
    
    % Remove all thresholded pixels.
    if remMode < 2
        
        % Dilate the thresholded pixels.
        if ~isempty(dilateDisk)
            bImg = imdilate(bImg, dilateDisk);
        end
        
        % Remove the thresholded pixels.
        img = single(img);
        img(bImg) = NaN;
        
        % Remove the worm.
    else
        
        % Find the worm.
        cc = bwconncomp(bImg);
        if ~isempty(cc.PixelIdxList)
            lengths = cellfun(@numel, cc.PixelIdxList);
            [~, ccMaxI] = max(lengths, [], 2);
            
            % Use the largest connected component.
            wormPixels = cc.PixelIdxList{ccMaxI};
            
            % Create a mask of the worm.
            mask = falseImg;
            mask(wormPixels) = true;
            
            % Dilate the worm.
            if ~isempty(dilateDisk)
                mask = imdilate(mask, dilateDisk);
            end
            
            % Remove the worm.
            img = single(img);
            img(mask) = NaN;
        end
    end
end
end

