classdef video_info < handle
    %
    %   Class:
    %   seg_worm.video_info
    %   
    
    properties
        fps
        file_name_info %Class: seg_worm.video.file_name_info
    end
    
    methods
        function obj = video_info(norm_folder)
           
            h = load(fullfile(norm_folder,'segNormInfo.mat'));
            obj.fps = h.fps;
            
            %obj.file_name_info = seg_worm.video.file_name_info(h.myAviInfo.url);
            
        end
    end
    
end

%{

Stored in segNormInfo.mat


myAviInfo
        approxFrameNum: -1
                   bpp: 24
                fourcc: 'mjpg'
                   fps: 25.8398
        frameTimeoutMS: 3000
                handle: 1380431973
                height: 480
    nHiddenFinalFrames: 3
             numFrames: 4639
                plugin: 'DirectShow'
         preciseFrames: 300
             timeStamp: 0
                  type: 'rgb'
                   url: [1x100 char]
                 width: 640


%}