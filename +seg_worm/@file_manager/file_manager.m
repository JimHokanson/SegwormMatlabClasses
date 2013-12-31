classdef file_manager < handle
    %
    %   Class:
    %   seg_worm.file_manager
    %
    %   This class is meant to handle different files that are used during
    %   the video parsing. These files are associated with the original
    %   video file. These files have not been shared publically by MRC.
    
    properties
       info_file  %xml file with experiment information
       log_file   %csv file with stage locations
       
       %NOTE: We might want to allow many different versions
       %of the diff file, based on the different ways it can be computed.
       diff_file  %mat file with video differentiation, this is computed 
       %after the video has been created and is not an "original" file
       
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

