classdef skeleton_indices < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.skeleton_indices
    %
    %   This was originally created for feature processing. I found a lot
    %   of off by 1 errors in the feature processing.
    %
    %   Used in: (list is not comprehensive)
    %   --------------------------------------------------------
    %   - posture bends
    %   - posture directions
    %
    %   NOTE: These are hardcoded for now. I didn't find much use in trying
    %   to make this dynamic based on some maximum value.
    %
    %   Typical Usage:
    %   --------------------------------------------------------
    %   SI = seg_worm.skeleton_indices;
        
    properties (Constant)
        HEAD_INDICES = 1:8;
        NECK_INDICES = 9:16;
        MID_INDICES  = 17:33;
        HIP_INDICES  = 34:41;
        TAIL_INDICES = 42:49;
        
        HEAD_TIP_INDICES  = 1:4;
        HEAD_BASE_INDICES = 5:8;
        TAIL_BASE_INDICES = 41:45;
        TAIL_TIP_INDICES  = 46:49;
        
% % %         %For getLocomotionBends
% % %         NOSE_INDICES = 4:-1:1 %Might remove this
% % %         NECK_INDICES = 8:-1:5 %Might remove this
        
        ALL_INDICES = 1:49   %Same as 1:end,
        
        BODY_INDICES = 9:41; %Indices excluding the head and tail
    end
    
    properties (Dependent)
       ALL_NORMAL_INDICES
    end
    
    properties
       ALL_NORMAL_NAMES = {'head' 'neck' 'midbody' 'hips' 'tail'};
    end
    
    methods
        function value = get.ALL_NORMAL_INDICES(obj)
           value = {obj.HEAD_INDICES obj.NECK_INDICES obj.MID_INDICES obj.HIP_INDICES obj.TAIL_INDICES};  
        end
    end
    
end

