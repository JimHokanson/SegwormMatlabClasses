classdef info < handle
    %
    %   Class:
    %   seg_worm.stage.info
    %
    %   Your one stop shop for experiment info parameters related to the
    %   stage!
    
    properties
        %steps/micron.
        microns_X
        microns_Y
        %steps/pixel
        pixels_X
        pixels_Y
        
        delay_time %Minimum time between stage movements or similarly, the
        %maximum time it takes for a stage movement to complete. If the
        %delay is too small the stage movements become chaotic.
    end
    
    methods
        function obj = info(exp_info_file)
            xml = xmlread(exp_info_file);
            
            % Read the steps/microns.
            obj.microns_X = str2double(xmlReadTag(xml, ...
                'configuration.info.stage.steps.equivalent.microns.x'));
            obj.microns_Y = str2double(xmlReadTag(xml, ...
                'configuration.info.stage.steps.equivalent.microns.y'));
            
            % Read the steps/pixels.
            obj.pixels_X = str2double(xmlReadTag(xml, ...
                'configuration.info.stage.steps.equivalent.pixels.x'));
            obj.pixels_Y = str2double(xmlReadTag(xml, ...
                'configuration.info.stage.steps.equivalent.pixels.y'));
            obj.delay_time = str2double(xmlReadTag(xml, ...
                'configuration.info.tracker.algorithm.delay'))/1000;
        end
    end
    
end

