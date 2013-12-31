classdef posture < handle
    %
    %   Class:
    %   seg_worm.features.posture
    
    properties
        bends
        amplitude
        wavelength
        trackLength
        eccentricity
        kinks
        coils
        directions
        skeleton
        eigenProjection  %[6 n_frames]
    end
    
    methods
        function obj = posture(nw,midbody_distance,FPS,p_opts)
                        
            %Bends
            %----------------------------------------------
            obj.getPostureBends(nw.angles);
            
            
            %Eccentricity & Orientation
            %---------------------------------------------------------
            worm_orientation = obj.getEccentricity(nw.contour_x, nw.contour_y, p_opts.N_ECCENTRICITY);
            
            
            %Amplitude, Wavelengths, TrackLength
            %--------------------------------------------------------
            obj.getAmplitudeAndWavelength(worm_orientation,nw.x,nw.y,nw.lengths);
            
            
            %Kinks (aka bend counts)
            %--------------------------------------------------------------
            obj.getWormKinks(nw.angles);
            
            
            %Coils
            %--------------------------------------------------------------------------
            
            %??? - Would we be able to identify coiled worms by integrating the angles?
            %cumsum(angles) - ensure cumsum never exceeds some value
            
            obj.getCoils(nw.frame_codes,midbody_distance,FPS);
            
            
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
    end
    
end

