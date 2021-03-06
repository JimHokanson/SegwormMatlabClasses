classdef vignette < handle
    %
    %   Class:
    %   seg_worm.vignette
    %
    %   Improvements:
    %   -------------
    %   1) I think the old code allowed for vignette construction
    %   from the frame data if the file was not present.
    %   2) We could create the object but keep track of a null case
    %   inside this object. This would allow the top level code
    %   to always 'apply' a vignette
    
    
    properties
        file_path
        mask_data %single, +/- values present
        
        %I think these are just from debugging
        img_scale = 1;
        img_data  %
    end
    methods
        function value = get.img_data(obj)
            if isempty(obj.img_data)
                obj.populateImageData();
            end
            value = obj.img_data;
        end
        function set.img_scale(obj,value)
           obj.img_scale = value;
           obj.populateImageData();
        end
    end
    
    methods (Static)
        function obj_or_empty = create(file_manager,video_reader)
            %
            %   vignette = seg_worm.vignette.create(file_manager,video_reader)
            %
            %   Inputs:
            %   -------
            %   file_manager : seg_worm.file_manager
            %   video_reader : seg_worm.videoReader
            %
            %   Outputs:
            %   --------
            %   vignette : 
            %
            %   NOTE: The static function allows us to return an empty object.
           v_file_path = file_manager.vignette_data_file;
           if exist(v_file_path,'file')
               obj_or_empty = seg_worm.vignette(v_file_path,video_reader);
           else
               %NOTE: I think this could automatically be created in the
               %old version
               obj_or_empty = [];
           end
        end
    end
    
    methods (Access = private)
        function obj = vignette(file_path,video_reader)
            %
            %   vignette = seg_worm.vignette.create(file_path,video_reader)
            %
            obj.file_path      = file_path;

            height = video_reader.height;
            width  = video_reader.width;

            %???? Why is the vignette signed????
            %
            %???? Why int32 to int8
            %
            %For some reason the data is written as a 4 byte value
            %when only occuping 1 byte
            temp = single(sl.io.fileRead(obj.file_path ,'int32=>int8','endian','b'));
            
            obj.mask_data = reshape(temp,width,height)';
        end
    end
    methods
        function fixed_data = apply(obj,frame_image_data)
            %
            %   Removes the vignette from an image frame
            %
            
            fixed_data = uint8(single(frame_image_data) - obj.mask_data);
            
            %Old code
            %uint8(single(fixed_images{iFrame}) - single(vImg));
        end
        
        function plot(obj)
            
            %mask
            subplot(1,2,1)
            imshow(obj.mask_data);
            
            % Show the vignette.
            subplot(1,2,2)
            imshow(obj.img_data);
            title(strrep(obj.file_path, '_', '\_')); % underscores confuse TeX
            
        end
    end
    methods (Hidden)
        function populateImageData(obj)
            
            % Color the positive values green and the negative values red.
            scale = obj.img_scale;
            
            vImg     = obj.mask_data/255;
            emptyImg = double(ones(size(vImg)));
            rImg = emptyImg;
            gImg = emptyImg;
            bImg = emptyImg;
            
            posImg       = vImg > 0;
            rImg(posImg) = 1 - vImg(posImg) * scale;
            
            bImg(posImg) = rImg(posImg);
            
            negImg       = vImg < 0;
            gImg(negImg) = 1 + vImg(negImg) * scale;
            bImg(negImg) = gImg(negImg);
            
            % Construct the RGB image.
            img(:,:,1) = rImg;
            img(:,:,2) = gImg;
            img(:,:,3) = bImg;
            
            obj.img_data = img;
            
        end
    end
    
end

