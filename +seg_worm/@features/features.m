classdef features < handle
    %
    %   Class:
    %   seg_worm.features
    %
    %   TODOS
    %   ===================================================================
    %   1) Implement info usage
    %   2) Implement saving and loading from struct
    %   3) Finish comparison function ...
    %   4) Implement options for processing and debugging ...
    
    properties
        morphology
        posture
        locomotion
        path
        info
    end
    
    methods
        function obj = features(nw,info,p_opts,d_opts)
            %
            %    seg_worm.features
            %
            %    Inputs
            %    ===================================================
            %    nw   : seg_worm.normalized_worm
            %    info : seg_worm.info
            %
            %    Optional Inputs
            %    ====================================================
            %    p_opts : seg_worm.features.processing_options
            %    d_opts : seg_worm.features.debug_options
            
            %TODO: Extract these from info
            FPS = 25.8398;
            VENTRAL_MODE = 0;  %??? I think this is set manually, but I'm not sure
            %where I should get this from at this point ...
            
            
            if ~exist('p_opts','var')
                p_opts = seg_worm.features.processing_options();
            end
            
            if ~exist('d_opts','var')
                d_opts = seg_worm.features.debug_options();
            end
            
            %processing_options - where we will specify how to process
            %things ...
            %debug_options - NYI
            
            obj.info = info;
            obj.morphology = seg_worm.features.morphology(nw,d_opts);
            obj.locomotion = seg_worm.features.locomotion(nw,FPS,VENTRAL_MODE,p_opts,d_opts);
            
            midbody_distance = abs(obj.locomotion.velocity.midbody.speed/FPS);
            
            obj.posture    = seg_worm.features.posture(nw,midbody_distance,FPS,p_opts,d_opts);
            
            obj.path = seg_worm.features.path(nw,FPS,VENTRAL_MODE);
            
        end
    end
    
end

