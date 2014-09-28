%function testWormVelocity()

addpath('C:\Users\mcurrie\Desktop\GitHub\SegwormMatlabClasses')

worm_path = 'C:\Backup\Dropbox\worm_data\video\testing_with_GUI\results\mec-4 (u253) off food x_2010_04_21__17_19_20__1_features.mat'
%'/Users/jameshokanson/Desktop/worm_data/segworm_data/features/T27E9.9 (ok2371)II on food R_2011_08_11__11_58_57___1___7_features.mat';


h = load(worm_path);

x = h.worm.posture.skeleton.x;
y = h.worm.posture.skeleton.y;

fps = 20;

velocity = seg_worm.features.locomotion.getWormVelocity(x, y, fps, 0);

midbody_speed = velocity.midbody.speed;

lengths = h.worm.morphology.length;

locomotion.motion = seg_worm.feature_helpers.locomotion.getWormMotionCodes(midbody_speed, lengths,fps);
