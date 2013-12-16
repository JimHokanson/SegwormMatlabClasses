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

fmp - feature mat path

seg_worm.feature_calculator.get_features_rewritten(NORM_PATH,fmp)

%}

%NOTE: I want seg_worm.normalized_worm to eventually
%be the input to this function. 

t_final = tic;

nw = seg_worm.normalized_worm.getObject(norm_folder);
%Class: seg_worm.normalized_worm

worm       = struct;


FPS = 25.8398;
VENTRAL_MODE = 0;  %??? I think this is set manually, but I'm not sure
%where I should get this from at this point ...
%
%   ventralMode:
%   0 - unknown
%   1 - clockwise
%   2 - anticlockwise
%
%   Used in:
%   seg_worm.feature_helpers.locomotion.getWormVelocity
%   seg_worm.feature_helpers.locomotion.getForaging
%

%Morphology - DONE
%--------------------------------------------------------------------------
worm.morphology = seg_worm.feature_calculator.getMorphologyFeatures(nw);


%Locomotion - DONE
%--------------------------------------------------------------------------
worm.locomotion = seg_worm.feature_calculator.getLocomotionFeatures(nw,FPS,VENTRAL_MODE);


%Posture
%--------------------------------------------------------------------------
midbody_distance = abs(worm.locomotion.velocity.midbody.speed/FPS);
worm.posture = seg_worm.feature_calculator.getPostureFeatures(nw,midbody_distance,FPS);


%Path
%--------------------------------------------------------------------------
worm.path = seg_worm.feature_calculator.getPathFeatures(nw,FPS,VENTRAL_MODE);

fprintf('Total Run Time %0.2g\n',toc(t_final));

%seg_worm.feature_calculator.verifyResult(worm,feature_mat_path);

