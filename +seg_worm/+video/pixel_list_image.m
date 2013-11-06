classdef pixel_list_image < handle
    %
    %   Class:
    %   seg_worm.video.pixel_list_image
    %
    %   Goal is to have a set of lists
    %   for an image which can be toggled.
    %
    %   Current Status:
    %   This currently helps with plotting results but I have yet to
    %   enable the turning off and on of layers by mouse-click
    %
    %
    
    properties
       default_color = [0 0 0]
       lists
       color
       names
       min_x = NaN
       max_x = NaN
       min_y = NaN
       max_y = NaN
    end
    
    methods
        function addList(obj,name,color,list_values)
           obj.lists = [obj.lists {list_values}];
           obj.color = [obj.color; color];
           obj.names = [obj.names {name}];
           
           obj.min_x = min(obj.min_x,min(list_values(:,2)));
           obj.min_y = min(obj.min_y,min(list_values(:,1)));
           obj.max_x = min(obj.max_x,max(list_values(:,2)));
           obj.max_y = min(obj.max_y,max(list_values(:,1)));
        end
        function plot(obj)
           %This is a work in progress ...
           final_image = zeros(obj.max_y,obj.max_x,3);
           for iList = 1:length(obj.lists)
              cur_list  = obj.lists{iList};
              cur_color = reshape(obj.color(iList,:),1,1,3);
              
              for iPixel = 1:size(cur_list,1)
                 I = cur_list(iPixel,1);
                 J = cur_list(iPixel,2); 
                 final_image(I,J,:) = cur_color;
              end
           end
           imshow(final_image);
        end
    end
    
end

