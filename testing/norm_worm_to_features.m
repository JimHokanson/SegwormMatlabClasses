%Normalized worm to feature file

%Merge blocked worm versions ...

NORM_PATH = 'F:\worm_data\segworm_data\video\testing_with_GUI\.data\mec-4 (u253) off food x_2010_04_21__17_19_20__1_seg\normalized';

%seg_worm.normalized_worm.createObjectFromFiles(NORM_PATH);

%nw = seg_worm.normalized_worm.getObject(NORM_PATH);

seg_worm.feature_calculator.get_features_rewritten(norm_folder)