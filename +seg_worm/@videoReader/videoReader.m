classdef videoReader
    %
    %   Class:
    %   seg_worm.videoReader
    %
    %   This is an interface class for the seg_worm code to
    %   a video reader. It may change as we change the video
    %   readers we use but it summarizes the things that need
    %   to be exposed in order for the seg_worm code to use
    %
    %   NOTE: The name might change to make it more obvious
    %   that this isn't a full fledged video reader
    
    properties
       r %sl.video.avi.reader
    end
    
    properties
       fps
       height
       width
       n_frames
    end
    
    methods
        function obj = videoReader(file_path,read_as_gray)
           %
           %
           %    
            
           %NOTE: sl is a reference to:
           %https://github.com/JimHokanson/matlab_standard_library
           obj.r      = sl.video.avi.reader(file_path);
           
           %NOTE: This is currently assuming that we are using
           %the mjpg codec where I can read only the luminance
           %values
           
           %This little bit allows me to speed up reading of the
           %the video
           if exist('read_as_gray','var') && read_as_gray
              obj.r.codec.read_as_gray = true;
           end
           
           r_local    = obj.r;
           obj.fps    = r_local.fps;
           obj.height = r_local.height;
           obj.width  = r_local.width;
           obj.n_frames = r_local.n_frames;
        end
        function [original_image,grayscale_image] = getFrame(obj,frame_number)
           %
           %
           %    
           original_image = obj.r.getFrame(frame_number);
           
           %Check for 
           if ismatrix(original_image)
              grayscale_image = original_image; 
           else
              %TODO: We might really want rgb2gray
              grayscale_image = original_image(:,:,1); 
           end
        end
        function close(obj)
           close(obj.r); 
        end
    end
    
end

