classdef file_manager < handle
    %
    %   Class:
    %   seg_worm.file_manager
    
    properties
       info_file  %xml file with experiment information
       log_file   %csv file with stage locations
       
       %NOTE: We might want to allow many different versions
       %of the diff file, based on different options ...
       diff_file  %mat file with video differentiation
       vignette_data_file %.dat file containing vignette info
    end
    
    methods
        function obj = file_manager(video_file_path)
           obj.vignette_data_file = strrep(video_file_path, '.avi', '.info.xml.vignette.dat');
           
           
           %NOTE: There are different options, we might want to allow
           %selecting between this different options ...
           obj.diff_file = strrep(video_file_path, '.avi', '_video_diff_file.mat');
           
           obj.log_file  = strrep(video_file_path, '.avi', '.log.csv');
           
           obj.info_file = strrep(video_file_path, '.avi', '.info.xml');
        end
    end
    
end

