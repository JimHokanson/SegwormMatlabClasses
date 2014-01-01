classdef locomotion < handle
    %
    %   Class:
    %   seg_worm.features.locomotion
    
    properties
        velocity
        %   .headTip
        %       .speed
        %       .direction
        %   .head
        %       .speed
        %       .direction
        %   .midbody
        %       .speed
        %       .direction
        %   .tail
        %       .speed
        %       .direction
        %   .tailTip
        %       .speed
        %       .direction
        motion
        bends
        turns
    end
    
    methods
        function obj = locomotion(nw,FPS,VENTRAL_MODE,p_opts,d_opts)
            
            %Velocity
            %--------------------------------------------------------------
            obj.getWormVelocity(nw.x,nw.y,FPS,VENTRAL_MODE,d_opts);
            
            
            %Motion
            %--------------------------------------------------------------
            midbody_speed     = obj.velocity.midbody.speed;
            obj.getWormMotionCodes(midbody_speed, nw.lengths, FPS);
            
            
            %Crawling - part of bends
            %--------------------------------------------------------------
            motion_mode = obj.motion.mode;
            is_paused   = motion_mode == 0;
            %SLOW - many ffts being computed
            obj.getLocomotionBends(nw.angles, is_paused, nw.is_segmented, FPS);
            
            
            %Foraging - part of bends
            %--------------------------------------------------------------
            obj.getForaging(nw.x,nw.y,nw.is_segmented,VENTRAL_MODE,FPS);
            
            
            %Turns
            %--------------------------------------------------------------
            midbody_distance  = abs(midbody_speed/FPS);
            is_stage_movement = nw.segmentation_status == 'm';
            
            obj.getOmegaAndUpsilonTurns(nw.angles,is_stage_movement,midbody_distance,nw.x,nw.y,FPS);
        end
    end
    
end

