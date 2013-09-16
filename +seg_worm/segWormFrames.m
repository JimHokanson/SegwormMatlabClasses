function [worms,orig_images,fixed_images] = segWormFrames(video_file_path, frames, verbose)
%SEGWORMFRAMES  Segment the worm in a set of video frames and organize
%   the information in a structure.
%
%   [worms,orig_images,fixed_images] = segWormFrames(videoFile, frames, verbose)
%
%   INPUTS
%   -----------------------------------------------------------
%   - At some point I lost 'num_erode' and 'num_dilate'
%   - I removed 'samples' and 'isInterp'
%   - See: seg_worm.worm.initialize
%
%   Inputs:
%       videoFile - the name of the video to segment
%       frames    - the frames of the video to segment
%                   Note: video frame indexing begins at 0
%                    use [] to allow all frames
%       verbose   - verbose mode shows the results in figures
%
%   Output:
%       worms - the worms' information organized in a cell array of structures
%               This structure contains 8 sub-structures,
%               6 sub-sub-structures, and 4 sub-sub-sub-structures:
%
%               * Video *
%               video = {frame}
%
%               * Contour *
%               contour = {pixels, touchI, inI, outI, angles, headI, tailI}
%
%               * Skeleton *
%               skeleton = {pixels, touchI, inI, outI, inOutI, angles,
%                           length, chainCodeLengths, widths}
%
%               * Head *
%               head = {bounds, pixels, area,
%                       cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%               head.bounds{contour.left (indices for [start end]),
%                           contour.right (indices for [start end]),
%                           skeleton indices for [start end]}
%
%               * Tail *
%               tail = {bounds, pixels, area,
%                       cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%               tail.bounds{contour.left (indices for [start end]),
%                           contour.right (indices for [start end]),
%                           skeleton indices for [start end]}
%
%               * Left Side (Counter Clockwise from the Head) *
%               left = {bounds, pixels, area,
%                       cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%               left.bounds{contour (indices for [start end]),
%                           skeleton (indices for [start end])}
%
%               * Right Side (Clockwise from the Head) *
%               right = {bounds, pixels, area,
%                        cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%               right.bounds{contour (indices for [start end]),
%                            skeleton (indices for [start end])}
%
%               * Orientation *
%               orientation = {head, vulva}
%               orientation.head = {isFlipped,
%                                   confidence.head, confidence.tail}
%               orientation.vulva = {isClockwiseFromHead,
%                                   confidence.vulva, confidence.nonVulva}
%
%       imgs  - a cell array of the requested frames converted to grayscale
%               and corrected for vignetting
%       oImgs - a cell array of the requested frames
%
%   See also WORM2STRUCT, SEGWORM, and SEGWORMVIDEO
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.



vr = videoReader(video_file_path,true);

if ~exist('frames','var') || isempty(frames)
   frames = 1:vr.n_frames;
end

if any(diff(frames) ~= 1)
    error('Non-consecutive frame requests not yet supported')
end

if ~exist('verbose','var') || isempty(verbose)
   verbose = true; 
end

file_manager = seg_worm.file_manager(video_file_path);

vignette     = seg_worm.vignette.create(file_manager,vr);
use_vignette = ~isempty(vignette);

% Pre-allocate memory.
orig_images  = cell(length(frames),1);
fixed_images = cell(length(frames),1);
worms        = cell(length(frames),1);

% Segment the video frames.
%--------------------------------------------------------------------------
n_frames    = length(frames);
worm_parsed = false(1,n_frames);
for iFrame = 1:n_frames
    cur_frame_number = frames(iFrame);
    
    %Retrieve frame data
    [oimg_temp,gimg_temp] = vr.getFrame(cur_frame_number);
    
    % Segment the worm and store its information.
    if isempty(oimg_temp)
%         warning('segWormFrames:NoFrame',...
%            'Cannot find frame %d, Perhaps it was dropped',...
%             cur_frame_number);
    else
        orig_images{iFrame} = oimg_temp;
        
        if use_vignette
            fixed_images{iFrame} = vignette.apply(gimg_temp);
        else
            fixed_images{iFrame} = gimg_temp;
        end
        
        temp_worm = seg_worm.worm(fixed_images{iFrame}, cur_frame_number, 1, verbose);
        
        worm_parsed(iFrame) = ~temp_worm.parse_error;
        worms{iFrame} = temp_worm;
        %NOTE: We might want to delete the original image which
        %is currently stored in the worm so that we can save memory
        
    end
end
% Clean up.
close(vr);

%JAH TODO: I'm at this point, need to implement code that is commented out
%below. This is all very confusing. I need to look back at what the overall
%strategy is for orienting the worm

for iFrame = 1:n_frames-1
   previous_worm  = worms{iFrame};
   cur_worm = worms{iFrame+1};
    
   if ~(previous_worm.parse_error || cur_worm.parse_error)
       %seg_worm.worm.orientation.orientToPreviousWorm
       cur_worm.orientation_handler.orientToRefWorm(previous_worm,verbose);
       
%        [hOrthoConfidence tOrthoConfidence ...
%  	hParaConfidence tParaConfidence ...
% 	hMagConfidence tMagConfidence] = ...
%     	headTailMovementConfidence(worms{iFrame}, next_worm);
       
       
       
   end
end



end

%This was in the loop following

%{

function helper__unused()
% Show the frame information.
    if ~isempty(worms{iFrame})
        hours   = floor(timestamp / 3600);
        minutes = floor((timestamp - hours * 60) / 60);
        seconds = (timestamp - hours * 3600 - minutes * 60);
        disp(['Worm at approximate frame = ' ...
            num2str(get(vr, 'approxFrameNum')) ...
            ', real frame = '  num2str(frame) ...
            ', timestamp = ' num2str(hours) ':' ...
            num2str(minutes, '%02.0f') ':' num2str(seconds, '%02.3f')]);
        
        % Compute the proximity and head/tail movement confidence.
        if next(vr)
            
            % Did we find the requested frame?
            timestamp = get(vr, 'timeStamp');
            frame = round(timestamp * fps);
            if frames(iFrame) + 1 ~= frame
                warning('segWormFrames:NoNextFrame', ['Frame ' ...
                    num2str(frames(iFrame) + 1) ' cannot be found. ' ...
                    'Therefore, the orientation and head/tail' ...
                    'movement confidences cannot be computed for its ' ...
                    'previous frame ' num2str(frames(iFrame))]);
                next_worm = [];
                
                % Get the next video frame and convert it to grayscale.
            else
                nextOImg = getframe(vr);
                nextImg = rgb2gray(nextOImg);
                nextImg = uint8(single(nextImg) - single(vImg));
                
                % Can the worm in the next frame be segmented?
                next_worm = segWorm(nextImg, frame, 1, verbose && ...
                    (iFrame < length(frames) && frames(iFrame + 1) == frame), ...
                    varargin{:});
                if isempty(next_worm)
                    warning('segWormFrames:NoNextWorm', ['Frame ' ...
                        num2str(frames(iFrame) + 1) ' cannot be segmented. ' ...
                        'Therefore, the orientation and head/tail' ...
                        'movement confidences cannot be computed for its ' ...
                        'previous frame ' num2str(frames(iFrame))]);
                    
                    % Orient the worm and compute the confidence.
                else
                    [next_worm confidence flippedConfidence] = ...
                        orientWorm(worms{iFrame}, next_worm, orientSamples);
                    [hOrthoConfidence tOrthoConfidence ...
                        hParaConfidence tParaConfidence ...
                        hMagConfidence tMagConfidence] = ...
                        headTailMovementConfidence(worms{iFrame}, next_worm);
                    
                    % Show the proximity and movement confidence.
                    disp(['Proximal orientation confidence:   ' 10 ...
                        '   Confidence = ' ...
                        num2str(confidence) 10 ...
                        '   Flipped confidence = ' ...
                        num2str(flippedConfidence)]);
                    disp(['Head/tail movement confidence: ' 10 ...
                        '   Head ortho-movement confidence = ' ...
                        num2str(hOrthoConfidence) 10 ...
                        '   tail ortho-movement confidence = ' ...
                        num2str(tOrthoConfidence) 10 ...
                        '   Head para-movement confidence = ' ...
                        num2str(hParaConfidence) 10 ...
                        '   Tail para-movement confidence = ' ...
                        num2str(tParaConfidence) 10 ...
                        '   Head movement magnitude confidence = ' ...
                        num2str(hMagConfidence) 10 ...
                        '   Tail movement magnitude confidence = ' ...
                        num2str(tMagConfidence)]);
                end
            end
        end
    end
end

%}