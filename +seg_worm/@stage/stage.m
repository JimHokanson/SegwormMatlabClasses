classdef stage < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.stage
    
    properties (Constant)
       MIN_FPS = 0.1;
       MAX_FPS = 100;
    end
    
    properties
       info           %seg_worm.stage.info

       fps      %the video frame rate (frames/second)
       frame_diffs  %the differences between subsequent video frames
       %See seg_worm.cv.video2Diff
        
       media_times
       locations
       
       %findStageMovement
       moves   = [0,0]
       origins = [0,0]
       
       %.readPixels2Microns()
       pixel_2_micron_scale = [-1, -1]
       rotation = 1
    end
    
    methods
        function obj = stage(file_manager)
           %
           %    INPUTS
           %    ----------------------------------
           %    file_manager : seg_worm.file_manager
           
           obj.info = seg_worm.stage.info(file_manager.info_file);
           
           [obj.media_times,obj.locations] = ...
               extractLogFileDetails(obj,file_manager.log_file);
           
           obj.readPixels2Microns();
           
           %TODO: If not present, we could run the diff code here ...
           h          = load(file_manager.diff_file);
           obj.fps        = h.fps;
           obj.frame_diffs = h.frameDiffs; 
           
           
           [frames,movesI,locations] = findStageMovement(obj);
           
           keyboard
           
        end
    end
    
end

