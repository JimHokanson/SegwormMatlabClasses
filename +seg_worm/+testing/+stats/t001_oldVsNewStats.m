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
    
    %TODO: Make this seg_worm.getHists for short
    hist_man_exp = seg_worm.stats.hist.manager(wormFiles(1:10));
    hist_man_ctl = seg_worm.stats.hist.manager(wormFiles(11:20));
    
    %TODO: Make this seg_worm.getStats for short
    stats_man = seg_worm.stats.manager(hist_man_exp,hist_man_ctl);
    
    hist_man.getStats();
    
    seg_worm.stats(hist_objs);
    
    %Old code for testing
    %----------------------------------------------------------------------
    for iFile = 1:10
        cur_file_path = wormFiles{iFile};
        [~,cur_name]  = fileparts(cur_file_path);
        cur_output_path = fullfile(hist_path,[cur_name '.mat']);
        seg_worm.w.stats.worm2histogram(cur_output_path, cur_file_path)
    end
    
    %Create a single hist file for controls ...
    tic
    cur_output_path = fullfile(hist_path,'control.mat');
    seg_worm.w.stats.worm2histogram(cur_output_path, wormFiles(11:20))
    toc
    
    
    %Merge all hist files into single file ...
    tic
    all_hist_files = cell(1,10);
    for iFile = 1:10
        cur_file_path = wormFiles{iFile};
        [~,cur_name]  = fileparts(cur_file_path);
        all_hist_files{iFile} = fullfile(hist_path,[cur_name '.mat']);

    end
    cur_output_path = fullfile(hist_path,'final_hist.mat');
    seg_worm.w.stats.addWormHistograms(cur_output_path, all_hist_files)
    toc
    
    tic
    test_path    = fullfile(hist_path,'final_hist.mat');
    control_path = fullfile(hist_path,'control.mat');
    merged_path  = fullfile(hist_path,'merged.mat');
    
    h1 = load(test_path);
    h2 = load(control_path);
    worm        = h1.worm;
    wormInfo    = h1.wormInfo;
    control     = h2.worm;
    controlInfo = h2.wormInfo;
    save(merged_path,'worm','wormInfo','control','controlInfo');
    
    %This doesn' work for some reason ...
    %Some bug with bin and edge specification
    seg_worm.w.stats.addWormHistograms(merged_path,test_path,control_path);
    toc
    

    stats_info_path = fullfile(hist_path,'stats_info.mat');
    %This approach is apparently incorrect ... :/
    %seg_worm.w.stats.worm2StatsInfo(stats_info_path, merged_path)
    
    seg_worm.w.stats.worm2StatsInfo(stats_info_path, test_path,[],[],control_path)
    
    seg_worm.stats.initObject(test_path, control_path)
    
%      controlData: [726x1 struct]
%      controlInfo: [1x10 struct]
%         dataInfo: [726x1 struct]
%     significance: [1x1 struct]
%         wormData: [726x1 struct]
%         wormInfo: [10x1 struct]
    
    
    
    
    
    
    