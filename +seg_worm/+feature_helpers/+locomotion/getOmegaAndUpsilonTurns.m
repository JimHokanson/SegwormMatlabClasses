function [omegas,upsilons] = getOmegaAndUpsilonTurns(bend_angles,is_stage_movement,midbody_distance,sx,sy,FPS)
%
%   JAH: I'm still working on this function ...
%
%   [omegas,upsilons] =
%   seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns(nw,tail_to_head_direction)
%
%   Old Name: code was in featureProcess.m and
%   omegaUpsilonDetectCurvature
%
%   Parent: 
%
%   See Also:
%   seg_worm.feature_calculator.getLocomotionFeatures   %Parent
%   seg_worm.feature_helpers.locomotion.getOmegaEvents
%   seg_worm.feature_helpers.locomotion.getTurnEventsFromSignedFrames
%   seg_worm.feature_helpers.locomotion.getUpsilonEvents
%   
%
%JAH: I'm at this point ...

%IMPORTANT: My events use 1 based indexing, their events use 0 based



INTER_DATA_SUM_NAME = 'interDistance';
DATA_SUM_NAME       = '';

n_frames = size(bend_angles,2);

SI = seg_worm.skeleton_indices;

%NOTE: For some reason the first and last few angles are NaN, so we use
%nanmean instead of mean
a.head_angles = nanmean(bend_angles(SI.FIRST_THIRD,:));
a.body_angles = nanmean(bend_angles(SI.SECOND_THIRD,:));
a.tail_angles = nanmean(bend_angles(SI.LAST_THIRD,:));

body_angles_for_ht_change = a.body_angles;

a.is_stage_movement = is_stage_movement;

n_head = sum(~isnan(a.head_angles));
n_body = sum(~isnan(a.body_angles));
n_tail = sum(~isnan(a.tail_angles));

%only proceed if there are at least two non-NaN value in each angle vector
if n_head < 2 || n_body < 2 || n_tail < 2
   %TODO: Create null event
   omegas   = seg_worm.feature.event.getNullStruct(FPS,DATA_SUM_NAME,INTER_DATA_SUM_NAME);
   upsilons = omegas;
   return 
end

%Interpolation
%--------------------------------------------------------------------------
a = h__interpolateAngles(a);

%Get frames for each turn type
%--------------------------------------------------------------------------
%This doesn't match was is written in the supplemental material ...
%Am I working off of old code??????
c = struct(...
    'headAngleStartConst',{20 -20 15 -15}, ...
    'tailAngleStartConst',{30 30  30  30}, ...
    'headAngleEndConst',  {40 40  30 30},  ...
    'tailAngleEndConst',  {20 -20 15 -15}, ...
    'bodyAngleConst'   ,  {20 -20 15 -15});

is_upsilon  = [false false true true];
sign_values = [1    -1     1    -1];

f.omegaFrames   = zeros(n_frames, 1);
f.upsilonFrames = zeros(n_frames, 1);

for iEntry = 1:4
    s = h__getConditionIndices(a,c(iEntry));
    f = h__populateFrames(a,s,f,is_upsilon(iEntry),sign_values(iEntry));
end

omegas   = seg_worm.feature_helpers.locomotion.getOmegaEvents(...
    f.omegaFrames,sx,sy,body_angles_for_ht_change,midbody_distance,FPS);
upsilons = seg_worm.feature_helpers.locomotion.getUpsilonEvents(...
    f.upsilonFrames,midbody_distance,FPS);

end

function fixed_x = h__interp_NaN(x)

fixed_x  = x;
nan_mask = isnan(x);

fixed_x(nan_mask) = interp1(find(~nan_mask),x(~nan_mask), find(nan_mask),'linear', 'extrap');

end

function a = h__interpolateAngles(a)

%TODO: Make the maximum gap explicit ...


%Get long NaN stretches ...
n = isnan(a.body_angles);
%This little bit finds runs of NaN values that are 10 samples or more
%0 -> A
%1 -> B
[long_nan_start_I, long_nan_end_I] = regexp( char(n+'A'), 'B{10,}', 'start', 'end' );

% interpolate arrays over NaN values (where there were stage
% movements, touching, or some other segmentation problem)
% ***This is of course only an approximate solution to the problem of
% not segmenting coiled shapes***
a.head_angles = h__interp_NaN(a.head_angles);
a.body_angles = h__interp_NaN(a.body_angles);
a.tail_angles = h__interp_NaN(a.tail_angles);

% return long NaN stretches back to NaN
for kk = 1:length(long_nan_start_I)
    a.bodyAngle(long_nan_start_I(kk):long_nan_end_I(kk)) = NaN;
end


end

function s = h__getConditionIndices(a,c)
%
%
%
%   This function implements a filter on the frames for the different
%   conditions that we are looking for in order to get a particular turn.
%   
%   It does not however provide any logic on their relative order. This is
%   done in a later function.

    %Determine comparison function
    %----------------------------------------------------------
    is_positive = c.headAngleStartConst > 0;

    if is_positive
        fh = @gt;
    else
        fh = @lt;
    end
    
    %start: when the head exceeds its angle but the tail does not
    %end  : when the tail exceeds its angle but the head does not
    
    s.startCond = fh(a.head_angles, c.headAngleStartConst) & abs(a.tail_angles) < c.tailAngleStartConst;
    s.startInds = find(diff(s.startCond) == 1) + 1; %add 1 for shift due to diff
    
    %NOTE: This is NaN check is a bit suspicious, as it implies that the
    %head and tail are parsed, but the body is not. The original code puts
    %NaN back in for long gaps in the body angle, so it is possible that
    %the body angle is NaN but the others are not.
    s.midCond   = fh(a.body_angles, c.bodyAngleConst) | isnan(a.bodyAngle);
    s.midStarts = find(diff(s.midCond) == 1) + 1; %add 1 for shift due to diff
    s.midEnds   = find(diff(s.midCond) == -1);
    
    s.endCond   = fh(a.tail_angles, c.tailAngleEndConst) & abs(a.head_angles) < c.headAngleEndConst;
    s.endInds   = find(diff(s.endCond) == -1);

end

function f = h__populateFrames(a,s,f,get_upsilon_flag,sign_value)
    
    %Algorithm:
    %-----------------------------------------------------------
    %- For the middle angle range, ensure one frame is valid and that
    %  the frame proceeding the start and following the end are valid
    %- Find start indices and end indices that bound this range
    %- For upsilons, exclude if they overlap with an omega bend ...
    
    %JAH NOTE: This type of searching is inefficient in Matlab since 
    %the data is already sorted. It could be improved ...
    
    for iMid = 1:length(s.midStarts)
        cur_mid_start_I = s.midStarts(iMid);
        cur_mid_end_I   = s.midEnds(find(s.midEnds > cur_mid_start_I, 1));
        
        if ~isempty(cur_mid_end_I)              && ...
            ~all(a.is_stage_movement(cur_mid_start_I:cur_mid_end_I)) && ...
            s.startCond(cur_mid_start_I - 1)    && ...
            s.endCond(cur_mid_end_I + 1)
 
            cur_start_I = s.startInds(find(s.startInds < cur_mid_start_I,   1, 'last'));
            cur_end_I   = s.endInds(find(s.endInds     > cur_mid_end_I,     1));

            if get_upsilon_flag
                %Don't populate upsilon if the data spans an omega
                if ~any(abs(f.omegaFrames(cur_start_I:cur_end_I)))
                    f.upsilonFrames(cur_start_I:cur_end_I) = sign_value;
                end
            else
                f.omegaFrames(cur_start_I:cur_end_I) = sign_value;
            end
        end
    end
end
