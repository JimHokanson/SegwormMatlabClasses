%function test_histogram_saving

%The code is really slow because a single field is loaded

% profile on

feature_root = 'F:\worm_data\segworm_data\features';

temp = sl.dir.getFilesInFolder(feature_root);

hist_output_path  = fullfile(feature_root,'results','hist.mat');
stats_output_path = fullfile(feature_root,'results','stats.mat');
final_stats_path  = fullfile(feature_root,'results','final_stats.mat');

%seg_worm.w.stats.worm2histogram(hist_output_path,temp.file_paths);

%seg_worm.w.stats.worm2stats(stats_output_path,hist_output_path);

%seg_worm.w.stats.worm2StatsInfo(final_stats_path,stats_output_path);


% profile off
% 
% profile viewer