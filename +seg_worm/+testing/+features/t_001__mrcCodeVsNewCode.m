function t_001__mrcCodeVsNewCode()
%
%
%   I want this function to test the old code vs my new code for computing
%   features.

%NORM_PATH = '/Users/jameshokanson/Dropbox/worm_data/video/testing_with_GUI/.data/mec-4 (u253) off food x_2010_04_21__17_19_20__1_seg/normalized'

base_path = 'F:\worm_data\segworm_data\video\testing_with_GUI';

NORM_PATH = fullfile(base_path,'.data\mec-4 (u253) off food x_2010_04_21__17_19_20__1_seg\normalized');
fmp = fullfile(base_path,'results\mec-4 (u253) off food x_2010_04_21__17_19_20__1_features.mat');

nw = seg_worm.normalized_worm.getObject(NORM_PATH);

info = seg_worm.info;

p_opts = seg_worm.features.processing_options;
d_opts = seg_worm.features.debug_options;

nf = seg_worm.features(nw,info,p_opts,d_opts);