function video_testing(use_new)
%Video testing
%----------------------
%
%   Current Status:
%
%   This is a test file for Jim which requires some videos that Ev sent to
%   me. 
%

if ispc
    base_path       = 'F:\worm_data\segworm_data\video';
else
    base_path       = '/Users/jameshokanson/Desktop/worm_data/segworm_data/video';
end
video_file_name = 'mec-4 (u253) off food x_2010_04_21__17_19_20__1.avi';
video_file_path = fullfile(base_path,video_file_name);

%vr = VideoReader(video_file_path);

%tic; temp = sl.media.mmread(video_file_path,1); toc;

% vr = sl.video.avi.reader(video_file_path);
% 
% for iFrame = 1:vr.n_frames
%    d = vr.getFrame(iFrame);
%    if ~isempty(d)
%       imshow(d)
%       drawnow
%    end
% end

% profile on


%I am currently examining the stage movement parsing ...
%failedFrames = seg_worm.saveWormFrames('wtf.mat', video_file_path);


%268 - worm is coiled - the code is not yet setup to handle this ...
FRAMES_PARSE = 1:250;

if use_new
tic
[worms, fixed_images, orig_images] = seg_worm.segWormFrames(video_file_path,FRAMES_PARSE,true);
toc
else
tic
[worms, fixed_images, orig_images] = segWormFrames(video_file_path,FRAMES_PARSE,true);
toc
end

% profile on
% tic
% seg_worm.cv.video2Diff(video_file_path);
% toc

% profile off
% profile viewer
