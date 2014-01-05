classdef posture < handle
    %
    %   Class:
    %   seg_worm.features.posture
    
    properties
        %Unless specified all sizes are: [1 n_frames]
        
        bends
        %   .head
        %       .mean
        %       .stdDev
        %   .neck
        %       .mean
        %       .stdDev
        %   .midbody
        %       .mean
        %       .stdDev
        %   .hips
        %       .mean     
        %       .stdDev  
        %   .tail
        %       .mean    
        %       .stdDev  
        amplitude 
        %   .max         
        %   .ratio       
        wavelength
        %   .primary     - 
        %   .wavelength  
        tracklength
        %   
        eccentricity
        
        kinks %(MF) (AKA Bend Counts) 
        
        coils %(EF)
        
        directions
        %   .tail2head
        %   .head
        %   .tail
        skeleton
        %   .x
        %   .y
        eigenProjection  %[6 n_frames]
    end
    
    methods
        function obj = posture(nw,midbody_distance,FPS,p_opts,d_opts)
                        
            %Bends
            %----------------------------------------------
            obj.getPostureBends(nw.angles,d_opts);
            
            
            %Eccentricity & Orientation
            %---------------------------------------------------------
            worm_orientation = obj.getEccentricity(nw.contour_x, nw.contour_y, p_opts.N_ECCENTRICITY);
            
            
            %Amplitude, Wavelengths, TrackLength
            %--------------------------------------------------------
            obj.getAmplitudeAndWavelength(worm_orientation,nw.x,nw.y,nw.lengths,d_opts);
            
            
            %Kinks (aka bend counts)
            %--------------------------------------------------------------
            obj.getWormKinks(nw.angles,p_opts);
            
            
            %Coils
            %--------------------------------------------------------------------------
            
            %??? - Would we be able to identify coiled worms by integrating the angles?
            %cumsum(angles) - ensure cumsum never exceeds some value
            
            obj.getCoils(nw.frame_codes,midbody_distance,FPS,d_opts);
            
            
            %Directions (AKA orientation)
            %--------------------------------------------------------------
            obj.getDirections(nw.x,nw.y);
            
            
            
            %Skeleton
            %--------------------------------------------------------------
            obj.skeleton.x = nw.x;
            obj.skeleton.y = nw.y;
            
            
            %EigenProjection
            %--------------------------------------------------------------
            obj.getEigenWorms(nw.x,nw.y,nw.eigen_worms,p_opts.N_EIGENWORMS_USE);
            
        end
        function plot(obj)
            
        end
    end
    
end

