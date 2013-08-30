classdef vignette < handle
    %
    %   Class:
    %   seg_worm.vignette
    
    properties
        img_data %single
    end
    
    methods
        function obj = vignette(vignette_file_path,height,width)
            
            %???? Why is the vignette signed????
            %
            %???? Why int32 to int8
            %
            %For some reason the data is written as a 4 byte value
            %when only occuping 1 byte
            temp = single(sl.io.fileRead(vignette_file_path,'int32=>int8','endian','b'));
            
            obj.img_data = reshape(temp,height,width)';
        end
        function fixed_data = apply(obj,frame_data)
            fixed_data = uint8(single(frame_data) - obj.img_data);
            
            %Old code
            %uint8(single(fixed_images{iFrame}) - single(vImg));
        end
    end
    
end

