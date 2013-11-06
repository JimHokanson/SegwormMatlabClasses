%Reproduction testing
%--------------------------------------

FEATURE_BASE_ROOT = 'F:\worm_data\segworm_data\features';

DATA_ROOT_1 = fullfile(FEATURE_BASE_ROOT,'acc-4','**','*.mat');
DATA_ROOT_2 = fullfile(FEATURE_BASE_ROOT,'gene_NA','**','*.mat');

d1 = seg_worm.fex.rdir(DATA_ROOT_1);
d2 = seg_worm.fex.rdir(DATA_ROOT_2);
