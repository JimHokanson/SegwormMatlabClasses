classdef stage < sl.obj.handle_light
    %
    %   Class:
    %   
    
    properties (Constant)
       MIN_FPS = 0.1;
       MAX_FPS = 100;
    end
    
    properties
       info_file_path %xml file containing experiment details
       log_file_path  %
       %findStageMovement
       moves   = [0,0]
       origins = [0,0]
       
       %.readPixels2Microns()
       pixel_2_micron_scale = [-1, -1]
       rotation = 1
    end
    
    methods
        function obj = stage(info_file_path,log_file_path)
           obj.info_file_path = info_file_path;
           obj.readPixels2Microns();
        end
    end
    
end

