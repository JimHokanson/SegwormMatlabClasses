function worm = get_features_rewritten(norm_folder)
%
%
%
%   seg_worm.feature_calculator.get_features_rewritten(norm_folder)


%{

TESTING CODE

NORM_PATH = 'F:\worm_data\segworm_data\video\testing_with_GUI\.data\mec-4 (u253) off food x_2010_04_21__17_19_20__1_seg\normalized';
NORM_PATH = '/Users/jameshokanson/Dropbox/worm_data/video/testing_with_GUI/.data/mec-4 (u253) off food x_2010_04_21__17_19_20__1_seg/normalized'
seg_worm.feature_calculator.get_features_rewritten(NORM_PATH)



%}

%NOTE: I want seg_worm.normalized_worm to eventually
%be the input to this function ...


%Load eigenworms
%load('masterEigenWorms_N2.mat');

nw = seg_worm.normalized_worm.getObject(norm_folder);
%Class: seg_worm.normalized_worm

worm       = struct;

%Morphology
%--------------------------------------------------------------------------
worm.morphology = seg_worm.feature_calculator.getMorphologyFeatures(nw);

%Locomotion
%--------------------------------------------------------------------------
%worm.locomotion = seg_worm.feature_calculator.getLocomotionFeatures(nw);

%NOTE: Also needs velocity as an input
%This is currently incomplete, I was waiting on some code
%from the MRC. They seem to have lost it. I'll probably need to rewrite
%it ...
%Posture
%--------------------------------------------------------------------------
worm.posture = seg_worm.feature_calculator.getPostureFeatures(nw);