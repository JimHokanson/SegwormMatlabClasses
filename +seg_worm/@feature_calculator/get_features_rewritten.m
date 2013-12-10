function worm = get_features_rewritten(norm_folder,feature_mat_path)
%
%
%
%   seg_worm.feature_calculator.get_features_rewritten(norm_folder)


%{
Awkward dependencies:
- coiling - seg_worm.feature_helpers.posture.getCoils - requires
worm.locomotion.velocity.midbody.speed


%}


%{

TESTING CODE
---------------------------------------------------------------------------

These files came from a video Ev gave to me. I then used the compiled worm
tracker GUI to generate the features. In the process of creating the
features I was able to generate intermediate files, of which the normalized
worm is one (or technically many, as it is divided into chunks of 500
frames per file)

NORM_PATH = 'F:\worm_data\segworm_data\video\testing_with_GUI\.data\mec-4 (u253) off food x_2010_04_21__17_19_20__1_seg\normalized';
NORM_PATH = '/Users/jameshokanson/Dropbox/worm_data/video/testing_with_GUI/.data/mec-4 (u253) off food x_2010_04_21__17_19_20__1_seg/normalized'
fmp = 'F:\worm_data\segworm_data\video\testing_with_GUI\results\mec-4 (u253) off food x_2010_04_21__17_19_20__1_features.mat';

seg_worm.feature_calculator.get_features_rewritten(NORM_PATH,fmp)

%}

%NOTE: I want seg_worm.normalized_worm to eventually
%be the input to this function. 

nw = seg_worm.normalized_worm.getObject(norm_folder);
%Class: seg_worm.normalized_worm

worm       = struct;

%Morphology - DONE
%--------------------------------------------------------------------------
worm.morphology = seg_worm.feature_calculator.getMorphologyFeatures(nw);

%Locomotion
%--------------------------------------------------------------------------
%worm.locomotion = seg_worm.feature_calculator.getLocomotionFeatures(nw);


%Posture
%--------------------------------------------------------------------------
worm.posture = seg_worm.feature_calculator.getPostureFeatures(nw);

seg_worm.feature_calculator.verifyResult(worm,feature_mat_path);

keyboard

%Path
%--------------------------------------------------------------------------
worm.path = seg_worm.feature_calculator.getPathFeatures(nw);