classdef processing_options < handle
    %
    %   Class:
    %   seg_worm.features.processing_options
    
    properties
       
    
    %Locomotion options
    %----------------------------------------------------------------------
    
        
    %Posture options
    %----------------------------------------------------------------------
    
    %seg_worm.features.posture.getEccentricity
    N_ECCENTRICITY = 50  %grid size for estimating eccentricity, this is the
    %max # of points that will fill the wide dimension
    
    %seg_worm.features.posture.getWormKinks
    KINK_LENGTH_THRESHOLD_PCT = 1/12; %This the fraction of the worm length
    %that a bend must be in order to be counted. Value is rounded to an 
    %integer sample. Threshold is inclusive.
    
    %seg_worm.features.posture.getEigenWorms
    N_EIGENWORMS_USE = 6 %# of eigenworms to extract. 
    %The maximum available number is 7.
    


    end
    
end

