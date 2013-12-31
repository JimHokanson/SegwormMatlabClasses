classdef processing_options < handle
    %
    %   Class:
    %   seg_worm.features.processing_options
    
    properties
       
    
    %Locomotion options
    %----------------------------------------------------------------------
    
        
    %Posture options
    %----------------------------------------------------------------------
    N_ECCENTRICITY = 50  %grid size for estimating eccentricity, this is the
    %max # of points that will fill the wide dimension
    
    N_EIGENWORMS_USE = 6 %# of eigenworms to extract. 
    %The maximum available number is 7.
    
    end
    
    methods
    end
    
end

