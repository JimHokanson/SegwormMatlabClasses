classdef videoReader
    %
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
           obj.r      = sl.video.avi.reader(file_path);
           
           %NOTE: This is currently assuming that we are using
           %the mjpg codec where I can read only the luminance
           %values
           
           if exist('read_as_gray','var') && read_as_gray
              obj.r.codec.read_as_gray = true;
           end
           
           r_local    = obj.r;
           obj.fps    = r_local.fps;
           obj.height = r_local.height;
           obj.width  = r_local.width;
           obj.n_frames = r_local.n_frames;
        end
        function [oimg,gimg] = getFrame(obj,frame_number)
           oimg = obj.r.getFrame(frame_number);
           
           if ismatrix(oimg)
              gimg = oimg; 
           else
              %TODO: We might really want rgb2gray
              gimg = oimg(:,:,1); 
           end
        end
        function close(obj)
           close(obj.r); 
        end
    end
    
end

