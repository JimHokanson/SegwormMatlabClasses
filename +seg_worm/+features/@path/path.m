classdef path < handle
    %
    %   Class:
    %   seg_worm.features.path
    %
    %   Nature Methods Description
    %   ===================================================================
    %   path.
    %   --------------
    %   The path features. 
    %
    %   The path is represented by its "range", "curvature", and the
    %   dwelling "duration" for various body parts. Individual experiment
    %   files also contain the "x" and "y" "coordinates" of the contour’s
    %   centroid.
    %   
    %   Moreover, the individual experiment files present the "duration" as
    %   an "arena" with a "height", "width", and the "min" and "max" values
    %   for the "x" and "y" axes of the arena. The arena can be transformed
    %   to a matrix using the given height and width.
    %
    %   The duration of the worm and body parts are represented as an array
    %   of "times" spent at the "indices" of the arena matrix
    
    
    properties (Hidden)
       nw 
    end
    
    properties
        range     %(MF) [1 x n_frames]
        
        duration  %(struct)
        %       arena: [1x1 struct]
        %           .height
        %           .width
        %           .min
        %              .x
        %              .y
        %           .max
        %              .x
        %              .y
        %        worm: 
        %           .indices
        %           .times   (SF)
        %        head:
        %           .indices
        %           .times   (SF)
        %     midbody:
        %           .indices
        %           .times   (SF)
        %        tail:
        %           .indices
        %           .times   (SF)
        %
        %   NOTE: Translation from indices to actual x & y locations is
        %   sort of simple but not just an indexing step ...

        coordinates
        %   .x - [1 x n_frames]
        %   .y - [1 x n_frames]
        
        curvature %[1 x n_frames]
    end
    
    methods
        function obj = path(nw,FPS,VENTRAL_MODE,d_opts)
            %
            %
            %
            
            obj.nw = nw;
            
            %Range
            %--------------------------------------------------------------
            obj.getRange(nw.contour_x,nw.contour_y);
            
            
            %Duration (aka Dwelling)
            %--------------------------------------------------------------
            obj.getDurationInfo(nw.x, nw.y, nw.widths, FPS,d_opts);
            
            
            %Coordinates
            %--------------------------------------------------------------
            obj.coordinates.x = mean(nw.contour_x);
            obj.coordinates.y = mean(nw.contour_y);
            
            
            %Curvature
            %--------------------------------------------------------------
            obj.wormPathCurvature(nw.x,nw.y,FPS,VENTRAL_MODE);
            
            
        end
        function plot(obj)
           keyboard 
           
           %Range
           %---------------------------------------------------------------
           figure(1)
           subplot(2,1,1)
           plot(obj.range,'Linewidth',3)
           title('Feature "range"')
           xlabel('Frame #','FontSize',18)
           ylabel('Range (microns)','FontSize',18)
           set(gca,'FontSize',18)
           subplot(2,1,2)
           x = obj.coordinates.x/10000;
           y = obj.coordinates.y/10000;
           scatter(nanmean(x),nanmean(y),100,'g','filled')
           hold on
           scatter(x,y,5,'r')
           
           valid_frames = find(~isnan(x));
           first_valid  = valid_frames(1);
           last_valid   = valid_frames(end);
           scatter(x(first_valid),y(first_valid),100,'k','filled')
           scatter(x(last_valid),y(last_valid),100,'b','filled')   
           hold off
           xlabel('X (cm)','FontSize',18)
           ylabel('Y (cm)','FontSize',18)
           set(gca,'FontSize',18)
           axis equal
           axis square
           legend({'Average of all centroids' 'Worm Centroids at each frame' '1st Frame' 'Last Frame' },'BestOutside')

           %Duration
           %---------------------------------------------------------------
           
        end
    end
    
end

