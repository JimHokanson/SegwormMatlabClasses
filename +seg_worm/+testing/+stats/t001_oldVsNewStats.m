function testing_code

    %seg_worm.stats.testing_code
 
    %seg_worm.w.stats.worm2histogram - create a hist file for a feature file
    
    
    %Newer code
    
    feature_path = '/Users/jameshokanson/Dropbox/worm_data/video/testing_with_GUI/results/mec-4 (u253) off food x_2010_04_21__17_19_20__1_features.mat';
    feature_path = 'F:\worm_data\segworm_data\video\testing_with_GUI\results\mec-4 (u253) off food x_2010_04_21__17_19_20__1_features.mat';
    
    
    
    base_path = 'F:\worm_data\segworm_data\features\gene_NA\allele_NA';
    hist_path = 'F:\worm_data\segworm_data\histograms';
    
    feature_files = sl.dir.rdir([base_path '\**\*.mat']);
    
    wormFiles = {feature_files.name};
    
    %This is somewhat arbitrary for now ...
    expt_files = wormFiles(1:10);
    ctl_files  = wormFiles(25:34);
    
    %TODO: Make this seg_worm.getHists for short
    hist_man_exp = seg_worm.stats.hist.manager(expt_files);
    hist_man_ctl = seg_worm.stats.hist.manager(ctl_files);
    
    %TODO: Make this seg_worm.getStats for short
    stats_man = seg_worm.stats.manager(hist_man_exp,hist_man_ctl);
    
    hist_man.getStats();
    
    seg_worm.stats(hist_objs);
    
    %TODO: This worked: write comparison script
    %stats_info_path_v3 = fullfile(hist_path,'stats_info_v3.mat');
    
    
    
    %Old code for testing
    %----------------------------------------------------------------------
    
    
    %APPROACH 1 ===========================================================
    hist_file = fullfile(hist_path,'all_worm_hists.mat');
    
    %NOTE: This approach might be really slow. Instead you could create
    %each hist individually then merge them later, but that takes more code
    %...
    %   NOTE: This doesn't seem to work
    worm2histogram(hist_file, expt_files, ctl_files);
    
    
    stats_info_path_v1 = fullfile(hist_path,'stats_info_v1.mat');
    %This approach is apparently incorrect ... :/
    %seg_worm.w.stats.worm2StatsInfo(stats_info_path, merged_path)
    
    %???? - pass in same file for both ?????
    %This doesn't work, error with concatenation ...
    %worm2StatsInfo(stats_info_path_v1, hist_file,[],[],hist_file)
    %
    %Same error
    %worm2StatsInfo(stats_info_path_v1, hist_file)
    
    %APPROACH 2  ==========================================================
    %Index exceeds matrix dimensions.
    
    expt_files = wormFiles(1:10);
    ctl_files  = wormFiles(25:34);
    
    expt_hist_file_old_v2 = fullfile(hist_path,'expt_hist_file_old_v2.mat');
    ctl_hist_file_old_v2  = fullfile(hist_path,'ctl_hist_file_old_v2.mat');
    
    worm2histogram(expt_hist_file_old_v2,expt_files);
    worm2histogram(ctl_hist_file_old_v2,ctl_files);
    
    stats_info_path_v2 = fullfile(hist_path,'stats_info_v2.mat');
    
    worm2StatsInfo(stats_info_path_v2, expt_hist_file_old_v2,[],[],ctl_hist_file_old_v2)
    
    %APPROACH 3 ===========================================================
    %This seems to be the approach actually used ...
    
    all_expt_hists_v3 = cell(1,10);
    for iFile = 1:10
        cur_file_path   = expt_files{iFile};
        [~,cur_name]    = fileparts(cur_file_path);
        cur_output_path = fullfile(hist_path,[cur_name '_v3.mat']);
        worm2histogram(cur_output_path, cur_file_path)
        all_expt_hists_v3{iFile} = cur_output_path;
    end
    
    all_ctl_hists_v3 = cell(1,10);
    for iFile = 1:10
        cur_file_path   = ctl_files{iFile};
        [~,cur_name]    = fileparts(cur_file_path);
        cur_output_path = fullfile(hist_path,[cur_name '_v3.mat']);
        worm2histogram(cur_output_path, cur_file_path)
        all_ctl_hists_v3{iFile} = cur_output_path;
    end
    
    expt_hist_merged = fullfile(hist_path,'expt_hist_file_old_merged_v3.mat');
    ctl_hist_merged  = fullfile(hist_path,'ctl_hist_file_old_merged_v3.mat');
    
    addWormHistograms(expt_hist_merged,all_expt_hists_v3);
    addWormHistograms(ctl_hist_merged,all_ctl_hists_v3);
    
    stats_info_path_v3 = fullfile(hist_path,'stats_info_v3.mat');
    
    worm2StatsInfo(stats_info_path_v3, expt_hist_merged,[],[],ctl_hist_merged)
    
    
    
    
    
    
    
    