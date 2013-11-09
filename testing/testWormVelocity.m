%function testWormVelocity()

worm_path = '/Users/jameshokanson/Desktop/worm_data/segworm_data/features/T27E9.9 (ok2371)II on food R_2011_08_11__11_58_57___1___7_features.mat';

h = load(worm_path);

x = h.worm.posture.skeleton.x;
y = h.worm.posture.skeleton.y;

fps = 20;

velocity = seg_worm.feature_helpers.locomotion.getWormVelocity(x, y, fps, 0);

midbody_speed = velocity.midbody.speed;

lengths = h.worm.morphology.length;

locomotion.motion = seg_worm.feature_helpers.locomotion.getWormMotionCodes(midbody_speed, lengths,fps);
