function getOmegaEvents(obj,omega_frames_from_angles,sx,sy,body_angles,midbody_distance,fps)
%
%
%
%   seg_worm.features.locomotion.getOmegaEvents
%
%   Inputs
%   =======================================================================
%   sx :
%   sy :
%   fps
%   body_angles : average bend angle of the middle third of the worm
%   midbody_distance :
%   omega_frames_from_angles : [1 x n_frames], each frame has the value 0,
%       1, or -1, 
%
%   Outputs
%   =======================================================================
%   omega_events : event structure 
%
%   Called By:
%   
%
%   See Also:
%   seg_worm.features.locomotion.getOmegaAndUpsilonTurns
%   seg_worm.features.locomotion.getTurnEventsFromSignedFrames

MIN_OMEGA_EVENT_LENGTH = round(fps/4);

body_angles_i = h__interp_NaN(body_angles,true); %_i - interpolated

omega_frames_from_th_change = h_getHeadTailDirectionChange(fps,sx,sy);

%Filter:
%This is to be consistent with the old code. We filter then merge, then
%filter again :/
omega_frames_from_th_change = h__filterAndSignFrames(...
    body_angles_i,omega_frames_from_th_change,MIN_OMEGA_EVENT_LENGTH);

is_omega_frame = omega_frames_from_angles | omega_frames_from_th_change;

%Refilter and sign
signed_omega_frames = h__filterAndSignFrames(body_angles_i,is_omega_frame,MIN_OMEGA_EVENT_LENGTH);

%Convert frames to events ...
obj.turns.omegas = obj.getTurnEventsFromSignedFrames(signed_omega_frames,midbody_distance,fps);

end

function is_omega_angle_change = h_getHeadTailDirectionChange(FPS,sx,sy)
%
%
%   NOTE: This change in direction of the head and tail indicates that
%   either a turn occurred OR that an error in the parsing occurred.
%   Basically we look for the angle from the head to the tail to all of a
%   sudden change by 180 degrees. 
%

MAX_FRAME_JUMP_FOR_ANGLE_DIFF = round(FPS/2);

%We compute a smoothed estimate of the angle change by using angles at
%indices that are +/- this value ...
HALF_WINDOW_SIZE = round(FPS/4);

%NOTE: It would be better to have this be based on time, not samples
MAX_INTERP_GAP_SIZE = 119;

%????!!!!?? - why is this a per frame value instead of an average angular
%velocity ????
PER_FRAME_DEGREE_CHANGE_CUTOFF = 3;


% Compute tail direction
%----------------------------------------------------
SI = seg_worm.skeleton_indices;

head_x = mean(sx(SI.HEAD_INDICES,:),1);
head_y = mean(sy(SI.HEAD_INDICES,:),1);
tail_x = mean(sx(SI.TAIL_INDICES,:),1);
tail_y = mean(sy(SI.TAIL_INDICES,:),1);

th_angle  = atan2(head_y - tail_y, head_x - tail_x)*(180/pi);

n_frames = length(th_angle);

%Changed angles to being relative to the previous frame
%--------------------------------------------------------------------------
%Compute the angle change between subsequent frames. If a frame is not
%valid, we'll use the last valid frame to define the difference, unless the
%gap is too large.

is_good_th_direction_value = ~isnan(th_angle);

lastAngle  = th_angle(1);
gapCounter = 0;

th_angle_diff_temp = NaN(size(th_angle));
for iFrame = 2:n_frames 
    if is_good_th_direction_value(iFrame)
        th_angle_diff_temp(iFrame) = th_angle(iFrame) - lastAngle;
        gapCounter = 0;
        lastAngle  = th_angle(iFrame);
    else
        gapCounter = gapCounter + 1;
    end
    if gapCounter > MAX_FRAME_JUMP_FOR_ANGLE_DIFF
        lastAngle = NaN;
    end
end

%???? - what does this really mean ??????
%I think this basically says, instead of looking for gaps in the original
%th_angle, we need to take into account how much time has passed between
%successive differences
%
%i.e. instead of doing a difference in angles between all valid frames, we
%only do a difference if the gap is short enough
positiveJumps = find(th_angle_diff_temp > 180);
negativeJumps = find(th_angle_diff_temp < -180);

%For example data, these are the indices I get ...
%P - 4625
%N - 3634, 4521 




%Fix the th_angles by unwrapping
%--------------------------------------------------------------------------
%NOTE: We are using the identified jumps from the fixed angles to unwrap
%the original angle vector
% subtract 2pi from remainging data after positive jumps
for j = 1:length(positiveJumps)
    th_angle(positiveJumps(j):end) = th_angle(positiveJumps(j):end) - 2*180;
end

% add 2pi to remaining data after negative jumps
for j = 1:length(negativeJumps)
    th_angle(negativeJumps(j):end) = th_angle(negativeJumps(j):end) + 2*180;
end



%Fix the th_angles through interpolation
%--------------------------------------------------------------------------
% get long NaN stretches
n = isnan(th_angle);
% save start and end indices for the stretches

gap_str = sprintf('B{%d,}',MAX_INTERP_GAP_SIZE+1); %Add 1 so that we allow the max gap
%but not anything greater
[start1, end1] = regexp( char(n+'A'), gap_str, 'start', 'end' );

% interpolate missing data
th_angle = h__interp_NaN(th_angle,false);

% return long NaN stretches back to NaN
for iEvent=1:length(start1)
    th_angle(start1(iEvent):end1(iEvent)) = NaN;
end

%Determine frames that might be omega events (we'll filter later based on
%length)
%--------------------------------------------------------------------------
% Compute angle difference
th_angle_diff = NaN(length(th_angle),1);

left_indices = (1:n_frames) - HALF_WINDOW_SIZE;
right_indics = (1:n_frames) + HALF_WINDOW_SIZE;

mask = left_indices > 1 & right_indics < n_frames;

th_angle_diff(mask) = th_angle(right_indics(mask)) - th_angle(left_indices(mask));

avg_angle_change_per_frame = abs(th_angle_diff/(HALF_WINDOW_SIZE*2));
is_omega_angle_change      = avg_angle_change_per_frame > PER_FRAME_DEGREE_CHANGE_CUTOFF;

end

function signed_omega_frames = h__filterAndSignFrames(body_angles_i,is_omega_frame,MIN_OMEGA_EVENT_LENGTH)

gap_str = sprintf('B{%d,}',MIN_OMEGA_EVENT_LENGTH);
[start1, end1] = regexp( char(is_omega_frame+'A')', gap_str, 'start', 'end');

signed_omega_frames = zeros(size(is_omega_frame));

%NOTE: Here we keep the long gaps instead of removing them

for iEvent = 1:length(start1)
    if mean(body_angles_i(start1(iEvent):end1(iEvent))) > 0
        signed_omega_frames(start1(iEvent):end1(iEvent)) = 1;
    else
        signed_omega_frames(start1(iEvent):end1(iEvent)) = -1;
    end
end

end


function fixed_x = h__interp_NaN(x,use_extrap)


fixed_x  = x;
nan_mask = isnan(x);

if use_extrap
    fixed_x(nan_mask) = interp1(find(~nan_mask),x(~nan_mask), find(nan_mask),'linear','extrap');
else
    fixed_x(nan_mask) = interp1(find(~nan_mask),x(~nan_mask), find(nan_mask),'linear');    
end

end
