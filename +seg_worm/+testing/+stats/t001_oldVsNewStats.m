function testing_code

    %seg_worm.stats.testing_code
 
    %seg_worm.w.stats.worm2histogram - create a hist file for a feature file
    
    
    %Newer code
    %----------------------------------------------------------------------
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
    
    stats_matrix_path_v3 = fullfile(hist_path,'stats_matrix_v3.mat');
    wormStats2Matrix(stats_matrix_path_v3, stats_info_path_v3)
    
    
    %Comparison code
    %----------------------------------------------------------------------
    %This needs to be moved somewhere else ...
    %
    %      control: [1x1 struct]
    %     dataInfo: [726x1 struct]
    %         worm: [1x1 struct]
    %
    %h.control.stats
    %        mean: [1x726 double]
    %      stdDev: [1x726 double]
    %     samples: [1x726 double]
    %     pNormal: [1x726 double]
    %     qNormal: [1x1 struct]
    %               strain: [1x726 double]
    %                  all: [1x726 double]  %This looks the same as strain
    %                  ....
    %
    %h.worm
    %      info: [1x1 struct]
    %     stats: [1x1 struct]
    %       sig: [1x1 struct]
    %   
    %h.worm.info
    %       strain: {'AQ2947'}
    %     genotype: {'C. elegans Wild Isolate, CGC N2 (Bristol, UK)'}
    %         gene: {[]}
    %       allele: {[]}
    %
    %h.worms.stats
    %
    %        mean: [1x726 double]
    %      stdDev: [1x726 double]
    %     samples: [1x726 double]
    %     pNormal: [1x726 double]
    %     qNormal: [1x1 struct]
    %      zScore: [1x726 double]
    %
    
    %????
    %- how does the order compare????

%
%                             m = morphology
%                             s = posture (shape)
%                             l = locomotion
%                             p = path
%
%                  type     = the feature's type, where:
%
%                             s = simple data
%                             m = motion data
%                             d = event summary data
%                             e = event data
%                             i = inter-event data
%
%                  subType  = the feature's sub-type, where:
%
%                             n = none
%                             f = forward motion data
%                             b = backward motion data
%                             p = paused data
%                             t = time data
%                             d = distance data
%                             h = frequency data (Hz)
%
%                  sign     = the feature's sign, where:
%
%                             s = signed data
%                             u = unsigned data
%                             a = the absolute value of the data
%                             p = the positive data
%                             n = the negative data
    
    
%     category: 'm'  - m,s,l,p
%         type: 'm'  - simple, motion, data, e vent, i nter-event
%      subType: 'n'  - none, forward, backward, paused, time, distance, h (freq)
%                       - this is for motion of the worm ...
%         sign: 'u'  - signed, unsigned, absolute, positive, negative
    
%     category: 'm'
%         type: 'm'
%      subType: 'f'
%         sign: 'u'

%     category: 'm'
%         type: 'm'
%      subType: 'p'
%         sign: 'u'

%     category: 'm'
%         type: 'm'
%      subType: 'b'
%         sign: 'u'


    %TODO: Need to reorganize my features to match order of their features
    %...

    
    stats_info_path_v3 = fullfile(hist_path,'stats_info_v3.mat');
    stats_matrix_path_v3 = fullfile(hist_path,'stats_matrix_v3.mat');
    
    h_si  = load(stats_info_path_v3);   %stats info
    h_sm  = load(stats_matrix_path_v3); %stats matrix
    h_stats = stats_man.stats;
    
    %h_si.dataInfo
    
    fields_all     = cell(1,726);
    categories_all = cell(1,726);
    type_all       = cell(1,726);
    sub_type_all   = cell(1,726);
    sign_all       = cell(1,726);
    for iObj = 1:726
       fields_all{iObj}     = h_si.dataInfo(iObj).field.histogram;
       categories_all{iObj} = h_si.dataInfo(iObj).category;
       type_all{iObj}       = h_si.dataInfo(iObj).type;
       sub_type_all{iObj}   = h_si.dataInfo(iObj).subType;
       sign_all{iObj}       = h_si.dataInfo(iObj).sign;
    end
    
        
    
    %add on .histogram.data.mean.all (or strip)
    
    %mean
    %mean.forward
    %mean.paused
    %mean.backward
    
    %.all
    %.abs
    %.pos
    %.neg
    
    %What about events
    
    field_names_new = cell
    
    
    'morphology.length.histogram.data.mean.all'
    'morphology.length.forward.histogram.data.mean.all'
    'morphology.length.paused.histogram.data.mean.all'
    'morphology.length.backward.histogram.data.mean.all'
    'morphology.width.head.histogram.data.mean.all'
    'morphology.width.head.forward.histogram.data.mean.all'
    'morphology.width.head.paused.histogram.data.mean.all'
    'morphology.width.head.backward.histogram.data.mean.all'
    'morphology.width.midbody.histogram.data.mean.all'
    'morphology.width.midbody.forward.histogram.data.mean.all'
    
    
    
    
    
    plot(h.worm.stats.zScore)
    hold all
    plot([stats_man.stats.z_score_experiment])
    hold off
    
    plot(h.worm.sig.pWValue)
    hold all
    plot([stats_man.stats.p_w])
    hold off
    
    %Similar but different ...
    plot(sort(h.worm.sig.pWValue))
    hold all
    plot(sort([stats_man.stats.p_w]))
    hold off
    
    %Manual comparison ...
    index = 1;
    h.worm.stats.pNormal(index)
    s(index).p_normal_experiment
    
    