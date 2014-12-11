function normalized_worm = parseWormFromVideo()
%
%   seg_worm.parseWormFromVideo()

BASE_PATH = 'C:\Users\RNEL\Dropbox\worm_data\video';
VERBOSE = true;
FRAMES = [];

video_name = 'mec-4 (u253) off food x_2010_04_21__17_19_20__1.avi';
video_path = fullfile(BASE_PATH,video_name);

%Main Old Code:
%- segWormVideo - NOTE: It is not setup to run. Any missing functions are
%either:
%1) Matlab toolboxes => e.g. image processing toolbox
%2) Can be found at:
%   @https://github.com/openworm/SegWorm
%3) Are 3rd party toolboxes ... :/

%Rough processing steps from the old code:
%1) progress_verbose.m - this is the main file that processes things:
%   2) segmentationMain.m - segment frames into skeletons and contours
%   3) headTail.m - head and tail detection
%   4) correctVulvalSide.m - correct the vulva side
%       OR
%      normWormProcess.m

[worms,orig_images,fixed_images] = seg_worm.segWormFrames(video_path, FRAMES, VERBOSE);

keyboard

end