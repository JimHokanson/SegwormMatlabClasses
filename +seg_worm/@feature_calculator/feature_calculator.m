classdef feature_calculator < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.feature_calculator
    %
    %   I might get rid of this class. It would be nice if features could
    %   calculate themselves.
    %
    %   
    %
    %
    %   I'm just throwing some processing code in here for now.
    
    properties
    end
    
    methods (Static)
        morphology = getMorphologyFeatures(nw)
        posture    = getPostureFeatures(nw,midbody_distance,FPS)
        locomotion = getLocomotionFeatures(nw,FPS,VENTRAL_MODE)
        path       = getPathFeatures(nw)
        get_features_rewritten(norm_folder,feature_mat_path)
        verifyResult(new_worm,feature_mat_path)
    end
    
end

