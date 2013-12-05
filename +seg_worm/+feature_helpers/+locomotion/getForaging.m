function foraging = getForaging(is_segmented_mask,sx,sy,ventralMode)
%
%
%   foraging = seg_worm.feature_helpers.locomotion.getForaging(is_segmented_mask,sx,sy,ventralMode)
%
%   Foraging
%   -----------------------
%   Foraging. Worm foraging is expressed as both an amplitude and an
%   angular speed (Supplementary Fig. 4g). Foraging is signed negatively
%   whenever it is oriented towards the ventral side. In other words, if
%   the nose is bent ventrally, the amplitude is signed negatively.
%   Similarly, if the nose is moving ventrally, the angular speed is signed
%   negatively. As a result, the amplitude and angular speed share the same
%   sign roughly only half the time. Foraging is an ambiguous term in
%   previous literature, encompassing both fine movements of the nose as
%   well as larger swings associated with the head. Empirically we have
%   observed that the nose movements are aperiodic while the head swings
%   have periodicity. Therefore, we measure the aperiodic nose movements
%   and term these foraging whereas the head swings are referred to as
%   measures of head crawling (described earlier in this section).
% 
%   Foraging movements can exceed 6Hz7 and, at 20-30fps, our video frame
%   rates are just high enough to resolve the fastest movements. By
%   contrast, the slowest foraging movements are simply a continuation of
%   the crawling wave and present similar bounds on their dynamics.
%   Therefore, we bound foraging between 1/30Hz (the lower bound used for
%   crawling) and 10Hz.
% 
%   To measure foraging, we split the head in two (skeleton points 1-4 and
%   5-8) and measure the angle between these sections. To do so, we measure
%   the mean of the angle between subsequent skeleton points along each
%   section, in the tail-to-head direction. The foraging angle is the
%   difference between the mean of the angles of both sections. In other
%   words, the foraging angle is simply the bend at the head.
% 
%   Missing frames are linearly interpolated, per each skeleton point, for
%   fragments up to 0.2 seconds long (4-6 frames at 20-30fps – twice the
%   upper foraging bound). When larger fragments are missing, foraging is
%   marked undefined. Segmentation of the head at very small time scales
%   can be noisy. Therefore, we smooth the foraging angles by convolving
%   with a Gaussian filter 1/5 of a second long (for similar reasons to
%   those mentioned in frame interpolation), with a width defined by the
%   Matlab “gausswin” function’s default ? of 2.5 and normalized such that
%   the filter integrates to 1.
% 
%   The foraging amplitude is defined as the largest foraging angle
%   measured, prior to crossing 0°. In other words, the largest nose bend
%   prior to returning to a straight, unbent position. Therefore, the
%   foraging amplitude time series follows a discrete, stair-step pattern.
%   The amplitude is signed negatively whenever the nose points towards the
%   worm’s ventral side. The foraging angular speed is measured as the
%   foraging angle difference between subsequent frames divided by the time
%   between these frames. To center the foraging angular speed at the frame
%   of interest and eliminate noise, each frame is assigned the mean of the
%   angular speed computed between the previous frame and itself and
%   between itself and the next frame. The angular speed is signed
%   negatively whenever its vector points towards the worm’s ventral side.
%
%   Improvements:
%   =======================================================================
%   - Change references to nose and neck to something more appropriate
%     NOTE: The neck is past the head, and the angle is within the head
%   - Change reference to nose_angles to something more appropriate

FPS = 20;

SI = seg_worm.skeleton_indices;
NOSE_I = SI.HEAD_TIP_INDICES(end:-1:1); %flip to maintain orientation for
%angles and consistency with old code ...
NECK_I = SI.HEAD_BASE_INDICES(end:-1:1);

MIN_NOSE_WINDOW = round(0.1*FPS);
MAX_NOSE_INTERP = 2*MIN_NOSE_WINDOW - 1;

nose_x = sx(NOSE_I,:);
nose_y = sy(NOSE_I,:);
neck_x = sx(NECK_I,:);
neck_y = sy(NECK_I,:);

interp_nose_mask = h__getNoseInterpolationIndices(is_segmented_mask,MAX_NOSE_INTERP);

dataI       = find(is_segmented_mask);
noseInterpI = find(interp_nose_mask);

nose_xi = h__interpData(nose_x,dataI,noseInterpI);
nose_yi = h__interpData(nose_y,dataI,noseInterpI);
neck_xi = h__interpData(neck_x,dataI,noseInterpI);
neck_yi = h__interpData(neck_y,dataI,noseInterpI);

noseBends = h__computeNoseBends(nose_xi,nose_yi,neck_xi,neck_yi);

[noseAmps,noseFreqs] = h__foragingData(noseBends, MIN_NOSE_WINDOW, FPS);

if ventralMode > 1
    noseAmps  = -noseAmps;
    noseFreqs = -noseFreqs;
end

foraging.amplitude  = noseAmps;
foraging.angleSpeed = noseFreqs;


end

function noseBends = h__computeNoseBends(nose_x,nose_y,neck_x,neck_y)

noseAngles = h__computeBendAngles(nose_x,nose_y);
neckAngles = h__computeBendAngles(neck_x,neck_y);

noseBends  = (noseAngles - neckAngles)';

wrap       = noseBends > pi;
noseBends(wrap) = noseBends(wrap) - 2 * pi;

wrap       = noseBends < -pi;
noseBends(wrap) = noseBends(wrap) + 2 * pi;

noseBends  = noseBends * 180 / pi; 

end

function angles = h__computeBendAngles(x,y)

    avg_diff_x = mean(diff(x,1,1));
    avg_diff_y = mean(diff(y,1,1));
    
    angles = atan2(avg_diff_y,avg_diff_x);

end


function x_new = h__interpData(x_old,dataI,noseInterpI)

x_new = x_old;

for i1 = 1:size(x_old,1)
   x_new(i1,noseInterpI) = interp1(dataI,x_old(i1,dataI),noseInterpI,'linear',NaN); 
end

end


function isInterpNoseData = h__getNoseInterpolationIndices(is_segmented_mask,maxNoseInterp)

%Interpolation ....
%==========================================================================
% Find the start and end indices for missing data chunks.
isNotData        = ~is_segmented_mask;
isInterpNoseData = isNotData;
diffIsNotData    = diff(isNotData);

%TODO: This needs to be fixed ...
startNotDataI    = find(diffIsNotData == 1);
endNotDataI      = find(diffIsNotData == -1);

% Don't interpolate missing data at the very start and end.
%--------------------------------------------------------------------------
if ~isempty(startNotDataI) && (isempty(endNotDataI) || startNotDataI(end) > endNotDataI(end))
    isInterpNoseData(startNotDataI(end):end) = false;
    startNotDataI(end) = [];
end

if ~isempty(endNotDataI) && (isempty(startNotDataI) || startNotDataI(1) > endNotDataI(1))
    isInterpNoseData(1:endNotDataI(1)) = false;
    endNotDataI(1) = [];
end

% Don't interpolate large missing chunks of data.

for i = 1:length(startNotDataI)
    if endNotDataI(i) - startNotDataI(i) > maxNoseInterp
        isInterpNoseData(startNotDataI(i):endNotDataI(i)) = false;
    end
end



end



%% Compute the foraging amplitude and angular speed.
function [amps,speeds] = h__foragingData(nose_bend_angle_d, minWinSize, fps)

% Clean up the signal with a gaussian filter.
if minWinSize > 0
    gaussFilter       = gausswin(2 * minWinSize + 1) / minWinSize;
    nose_bend_angle_d = conv(nose_bend_angle_d, gaussFilter, 'same');
    
    nose_bend_angle_d(1:minWinSize) = NaN;
    nose_bend_angle_d((end - minWinSize + 1):end) = NaN;
end


%NOTE: This code doesn't handle unsegmented data, so the foraging magnitude
%could be inaccurate in many places ...

%JAH NOTE: I left this code as is for now, I could improve it but
%it isn't that slow
%--------------------------------------------------------------------------
% Compute the amplitudes between zero crossings.
tic
dataSign = sign(nose_bend_angle_d);
amps     = NaN(1,length(nose_bend_angle_d));
numAmps  = 0;
for i = 1:(length(nose_bend_angle_d) - 1)
    
    % Compute the amplitude for the region.
    % Note: data at the zero crossing has NaN (unknown) amplitude.
    if dataSign(i) ~= dataSign(i + 1);
        if dataSign(i) > 0
            amps((i - numAmps):i) = max(nose_bend_angle_d((i - numAmps):i));
        elseif dataSign(i) < 0
            amps((i - numAmps):i) = min(nose_bend_angle_d((i - numAmps):i));
        end
        
        % Reset the count.
        numAmps = 0;
        
    % Advance.
    else
        numAmps = numAmps + 1;
    end
end

% Compute the amplitude for the end region.
% Note: data at the zero crossing has NaN (unknown) amplitude.
if dataSign(end) > 0
    amps((end - numAmps):end) = max(nose_bend_angle_d((end - numAmps):end));
elseif dataSign(end) < 0
    amps((end - numAmps):end) = min(nose_bend_angle_d((end - numAmps):end));
end
toc
%--------------------------------------------------------------------------


% Compute the speed centered between the back and front foraging movements.
%
%  1     2    3
%    d1    d2
%        x        assign to x, avg of d1 and d2
%
dData  = diff(nose_bend_angle_d) * fps;
speeds = NaN(size(amps));
speeds(2:end-1) = (dData(1:(end - 1)) + dData(2:end)) / 2;

amps(isnan(speeds)) = NaN;

end