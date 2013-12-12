function [omegas,upsilons] = getOmegaAndUpsilonTurns(bend_angles,is_stage_movement,midbody_distance,sx,sy,FPS)
%
%   seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns
%
%   Old Name: 
%   - featureProcess.m
%   - omegaUpsilonDetectCurvature.m
%
%   See Also:
%   seg_worm.feature_calculator.getLocomotionFeatures   %Parent
%   seg_worm.feature_helpers.locomotion.getOmegaEvents
%   seg_worm.feature_helpers.locomotion.getTurnEventsFromSignedFrames
%   seg_worm.feature_helpers.locomotion.getUpsilonEvents
%   
%   Nature Methods Description
%   =======================================================================
%   Turns.
%   ----------------------------------
%   Omega and upsilon turn events are computed similarly to a previously
%   described method9 but using skeleton bends instead of a single
%   head-midbody-tail angle. Omega and upsilon turns are signed negatively
%   whenever the worm’s ventral side is sheltered within the concavity of
%   its midbody bend.
%
%   The worm bends (described in the section on “Posture”) are used to find
%   a contiguous sequence of frames (interruptible by coiling and other
%   segmentation failures) wherein a large bend travels from the worm’s
%   head, through its midbody, to its tail. The worm’s body is separated
%   into three equal parts from its head to its tail. The mean
%   supplementary angle is measured along each third. For omega turns, this
%   angle must initially exceed 30° at the first but not the last third of
%   the body (the head but not the tail). The middle third must then exceed
%   30°. And finally, the last but not the first third of the body must
%   exceed 30° (the tail but not the head). This sequence of a 30° mean
%   supplementary angle, passing continuously along the worm from head to
%   tail, is labeled an omega turn event. Upsilon turns are computed nearly
%   identically but they capture all events that escaped being labeled
%   omega turns, wherein the mean supplementary angle exceeded 15° on one
%   side of the worm (the first or last third of the body) while not
%   exceeding 30° on the opposite end.

%IMPORTANT: My events use 1 based indexing, their events use 0 based

MAX_INTERPOLATION_GAP_ALLOWED = 9;
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
   omegas   = seg_worm.feature.event.getNullStruct(FPS,DATA_SUM_NAME,INTER_DATA_SUM_NAME);
   upsilons = omegas;
   return 
end

%Interpolation
%--------------------------------------------------------------------------
a = h__interpolateAngles(a,MAX_INTERPOLATION_GAP_ALLOWED);

%Get frames for each turn type
%--------------------------------------------------------------------------
%This doesn't match was is written in the supplemental material ...
%Am I working off of old code??????
c = struct(...
    'head_angle_start_const',{20 -20 15 -15}, ...
    'tail_angle_start_const',{30  30 30  30}, ...
    'head_angle_end_const',  {40  40 30  30}, ...
    'tail_angle_end_const',  {20 -20 15 -15}, ...
    'body_angle_const'   ,   {20 -20 15 -15});

%NOTE: We need to run omegas first (false values) since upsilons are more
%inclusive, but can not occur if an omega event occurs
is_upsilon  = [false false true true];

%NOTE: We assign different values based on the sign of the angles
values_to_assign = [1    -1     1    -1];

f.omegaFrames   = zeros(n_frames, 1);
f.upsilonFrames = zeros(n_frames, 1);

for iEntry = 1:4
    s = h__getConditionIndices(a,c(iEntry));
    f = h__populateFrames(a,s,f,is_upsilon(iEntry),values_to_assign(iEntry));
end

%Calculate the events from the frame values
%--------------------------------------------------------------------------
omegas   = seg_worm.feature_helpers.locomotion.getOmegaEvents(...
    f.omegaFrames,sx,sy,body_angles_for_ht_change,midbody_distance,FPS);
upsilons = seg_worm.feature_helpers.locomotion.getUpsilonEvents(...
    f.upsilonFrames,midbody_distance,FPS);

end

function fixed_x = h__interp_NaN(x)
%
%   TODO: Incorporate into 
%   seg_worm.feature_helpers.interpolateNanData
%

fixed_x  = x;
nan_mask = isnan(x);

fixed_x(nan_mask) = interp1(find(~nan_mask),x(~nan_mask), find(nan_mask),'linear', 'extrap');

end

function a = h__interpolateAngles(a,MAX_INTERPOLATION_GAP_ALLOWED)
%
%
%   Inputs
%   =======================================================================
%   a.head_angles
%   a.body_angles
%   a.tail_angles
%
%   Outputs
%   =======================================================================

%   TODO: Incorporate into 
%   seg_worm.feature_helpers.interpolateNanData


%Get long NaN stretches ...
n = isnan(a.body_angles);
%This little bit finds runs of NaN values that are 10 samples or more
%0 -> A
%1 -> B

str = sprintf('B{%d,}',MAX_INTERPOLATION_GAP_ALLOWED+1);

[long_nan_start_I, long_nan_end_I] = regexp( char(n+'A'), str, 'start', 'end' );

% interpolate arrays over NaN values (where there were stage
% movements, touching, or some other segmentation problem)
% ***This is of course only an approximate solution to the problem of
% not segmenting coiled shapes***
a.head_angles = h__interp_NaN(a.head_angles);
a.body_angles = h__interp_NaN(a.body_angles);
a.tail_angles = h__interp_NaN(a.tail_angles);

% return long NaN stretches back to NaN- only for the body angles ...
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
%   It does not however provide any logic on their relative order, i.e.
%   that one condition occurs before another. This is done in a later
%   function, h__populateFrames.

    %Determine comparison function
    %----------------------------------------------------------
    is_positive = c.head_angle_start_const > 0;

    if is_positive
        fh = @gt;
    else
        fh = @lt;
    end
    
    %start: when the head exceeds its angle but the tail does not
    %end  : when the tail exceeds its angle but the head does not
    
    s.startCond = fh(a.head_angles, c.head_angle_start_const) & abs(a.tail_angles) < c.tail_angle_start_const;
    s.startInds = find(diff(s.startCond) == 1) + 1; %add 1 for shift due to diff
    
    %NOTE: This is NaN check is a bit suspicious, as it implies that the
    %head and tail are parsed, but the body is not. The original code puts
    %NaN back in for long gaps in the body angle, so it is possible that
    %the body angle is NaN but the others are not.
    s.midCond   = fh(a.body_angles, c.body_angle_const) | isnan(a.bodyAngle);
    s.midStarts = find(diff(s.midCond) == 1) + 1; %add 1 for shift due to diff
    s.midEnds   = find(diff(s.midCond) == -1);
    
    s.endCond   = fh(a.tail_angles, c.tail_angle_end_const) & abs(a.head_angles) < c.head_angle_end_const;
    s.endInds   = find(diff(s.endCond) == -1);

end

function f = h__populateFrames(a,s,f,get_upsilon_flag,value_to_assign)
%
%
%   Inputs
%   =======================================================================
%    a: (structure)
%           head_angles: [1x4642 double]
%           body_angles: [1x4642 double]
%           tail_angles: [1x4642 double]
%     is_stage_movement: [1x4642 logical]
%             bodyAngle: [1x4642 double]
%    s: (structure)
%     startCond: [1x4642 logical]
%     startInds: [1x81 double]
%       midCond: [1x4642 logical]
%     midStarts: [268 649 881 996 1101 1148 1202 1963 3190 3241 4144 4189 4246 4346 4390 4457 4572 4626]
%       midEnds: [301 657 925 1009 1103 1158 1209 1964 3196 3266 4148 4200 4258 4350 4399 4461 4579]
%       endCond: [1x4642 logical]
%       endInds: [1x47 double]
%   f: (structure)
%       omegaFrames: [4642x1 double]
%     upsilonFrames: [4642x1 double]
%   get_upsilon_flag : toggled based on whether or not we are getting
%               upsilon events or omega events
%   sign_value : 
%
%   Outputs
%   =======================================================================
%
    
    %Algorithm:
    %-----------------------------------------------------------
    %- For the middle angle range, ensure one frame is valid and that
    %  the frame proceeding the start and following the end are valid
    %- Find start indices and end indices that bound this range
    %- For upsilons, exclude if they overlap with an omega bend ...
    

    
    for iMid = 1:length(s.midStarts)
        cur_mid_start_I = s.midStarts(iMid);
        
        %JAH NOTE: This type of searching is inefficient in Matlab since 
        %the data is already sorted. It could be improved ...
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
                    f.upsilonFrames(cur_start_I:cur_end_I) = value_to_assign;
                end
            else
                f.omegaFrames(cur_start_I:cur_end_I) = value_to_assign;
            end
        end
    end
end
