%testing

%1) Highest level processor
%2) Image level processor
%   - segWorm
%   

base_path = 'F:\worm_data\segworm_data';


video_file_name = 'mec-4 (u253) off food x_2010_04_21__17_19_20__1.avi';
video_file_path = fullfile(base_path,video_file_name);
worm_file_path  = sl.dir.changeFileExtension(video_file_path,'.mat');
tic
failedFrames    = saveWormFrames(worm_file_path, video_file_path);
toc