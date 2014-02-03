function t001_oldVsNewStats
%
%
%   Calling this code
%
%
%   NOTES: I am using the original code from MRC, not my modified versions
%   which have been cleaned up
%
%   BUGS FOUND:
%   ---------------------------------------
%   1) removePartialEvents - end match is wrong
%
%   Frequency - my defaults might be wrong - no events should be 0, not NaN
%   
%   This should be clearly stated in the events spec, rather than in the
%   code ...
%
%   MY ERRORS - that I haven't fixed yet
%   ================================================
%   1) I am not signing certain events properly
%   
%   - signing looks like it happens in worm2histogram() in event2histograms()
%
%   2) I need to add on that certain events have defaults of zero if the
%   event never occurred ...
%
%   - code might be in addWormHistograms() in addEventData()
%
%   HUH
%   -------------------------------------------------
%   1) In the Nature Methods paper, crawling features while the body motion
%   is paused are omitted. (Specifically amplitude and frequency for head,
%   midbody, and tail)
%
%   Feature Count
%   --------------------------------------------------
%   726 - old files
%   708 - new files (events need to be signed)
%   702 - nature methods (removal of paused body motion, 6 categories, 4
%   features each paused all, paused abs., paused neg, paused pos, 24
%   total)
%

    %seg_worm.stats.testing_code
    %seg_worm.w.stats.worm2histogram - create a hist file for a feature file
    
    %Current summary: The newer code needs to be run manually.
    
    
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



    %======================================================================
    %======================================================================
    %                     Comparison Code
    %======================================================================
    %======================================================================

    stats_info_path_v3 = fullfile(hist_path,'stats_info_v3.mat');
    stats_matrix_path_v3 = fullfile(hist_path,'stats_matrix_v3.mat');
    
    h_si    = load(stats_info_path_v3);   %stats info - old objects
    h_sm    = load(stats_matrix_path_v3); %stats matrix - old objects
    h_stats = stats_man.stats; %new objects
    
    %h_si.dataInfo
    
    h__generateOldStats
    
    [index_new__value_old,field_names_new,fields_old] = h_getOrder(h_si,h_sm,h_stats);

%==========================================================================    
%                   Comparison Code    
%==========================================================================
   
field_names_new_ordered(index_new__value_old) = field_names_new;
h_new(index_new__value_old) = h_stats;

clear h_exp_new
clear h_ctl_new
h_exp_new(index_new__value_old) = hist_man_exp.hists;
h_ctl_new(index_new__value_old) = hist_man_ctl.hists;

keep_mask    = ~cellfun('isempty',{h_new.name});
kept_indices = find(keep_mask);

%We don't have the mean or std data here
%mean control
%mean_control_old = h_sm.control.stats.mean(keep_mask);
%mean_control_new = [

%zscore_control - 0 by definition

%zscore_experiment
zscore_experiment_old = h_sm.worm.stats.zScore(keep_mask);
zscore_experiment_new = [h_new(keep_mask).z_score_experiment];

d_zscore_experiment = zscore_experiment_old - zscore_experiment_new;

I = find(abs(d_zscore_experiment) > 0.00001);


%Discrepancy: - AND MORE INDICES - some are NaN
%----------------------------------------------------------
%h_ctl_new(217).mean
%h_sm.control.stats.mean(217)

%Indices new
%677   681   685   691   695   699
%Indices old
%677   678   679   682   683   684

%     {h_ctl_new(kept_indices(I)).field}'
%
%
%   ???? Why is the order changing ?????
%
%     'locomotion.turns.omegas.frames.time'
%     'locomotion.turns.omegas.frames.interTime'
%     'locomotion.turns.omegas.frames.interDistance'
%     'locomotion.turns.upsilons.frames.time'
%     'locomotion.turns.upsilons.frames.interTime'
%     'locomotion.turns.upsilons.frames.interDistance'


%NOTES: 677 - all old are negative except one from experiment data, wtf?
%
%   678 - completely different scales - is resolution or something wrong?
%   679 - many from old are NaN (in fact almost all)
%
%   These might not be aligned correctly ...
%   NOPE: These are not aligned correctly ... WHY NOT?????

IDX = 682;
IDX_NEW = kept_indices(IDX);
fprintf(2,'Control Data\n');
sort(h_ctl_new(IDX_NEW).mean_per_video)
sort(h_si.controlData(IDX).dataMeans)
fprintf(2,'Experiment Data\n');
sort(h_exp_new(IDX_NEW).mean_per_video)
sort(h_si.wormData(IDX).dataMeans)


%Names:
%---------------------
h_si.dataInfo(IDX).field
h_ctl_new(IDX_NEW).field
     
%worm2histogram('wasfdsadf',ctl_files{5})

%motion.forward.frames.time
%motion.forward.frames.distance

%Let's look at the calculated statistics
%----------------------------------------------------------
pt_old  = h_sm.worm.sig.pTValue(keep_mask);
pt_new  = [h_new(keep_mask).p_t];
pt_diff = pt_old - pt_new;

%h_sm.worm.stats.zScore(I(1))
%h_new(1)

%h_sm.worm.stats.mean(217)
%h_sm.worm.stats.stdDev(217)
%h_sm.control.stats.mean(217)
%h_sm.control.stats.stdDev(217)


%Temporary Code
%==========================================================================
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
end

function [index_new__value_old,field_names_new,fields_old] = h_getOrder(h_si,h_sm,h_stats)


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


    fields_old     = cell(1,726);
    categories_all = cell(1,726);
    type_all       = cell(1,726);
    sub_type_all   = cell(1,726);
    sign_all       = cell(1,726);
    for iObj = 1:726
       fields_old{iObj}     = h_si.dataInfo(iObj).field.histogram;
       categories_all{iObj} = h_si.dataInfo(iObj).category;
       type_all{iObj}       = h_si.dataInfo(iObj).type;
       sub_type_all{iObj}   = h_si.dataInfo(iObj).subType;
       sign_all{iObj}       = h_si.dataInfo(iObj).sign;
    end
        
    field_names_new = cell(1,708);
    for iObj = 1:708       
       cur_stat = h_stats(iObj);
       cur_field_name = cur_stat.field;

       if strcmp(cur_stat.hist_type,'event')
           
           %Different rules for frames and no frames
           
           I = strfind(cur_field_name,'frames');
           
           if ~isempty(I)
           
           cur_field_name = regexprep(cur_field_name,'\.frames\.','\.');
           pre_str  = '';
           post_str = '.all';
           mid_str = '.histogram.data.mean';
           else
                pre_str  = '';
                mid_str = '.data';
                post_str = '';
           end
       elseif strcmp(cur_stat.hist_type,'simple')
           cur_field_name = regexprep(cur_field_name,'\.times','');
           pre_str  = '';
           post_str = '.all';
           mid_str = '.histogram.data.mean';
       else
           switch cur_stat.motion_type
               case 'all'
                   pre_str = '';
               otherwise
                   %forward, paused, backward
                   pre_str = ['.' cur_stat.motion_type];
           end

           %Indexing add on :/
           if sl.str.contains(cur_field_name,'eigenProjection')
              %We need to know the index ... 
              index_number = cur_stat.short_name(end); %Hack for getting #
              pre_str = ['(' index_number ')' pre_str]; %#ok<AGROW>
           end

           switch cur_stat.data_type
               case 'all'
                   post_str = '.all';
               case 'absolute'
                   post_str = '.abs';
               case 'negative'
                   post_str = '.neg';
               case 'positive'
                   post_str = '.pos';
           end
           mid_str = '.histogram.data.mean';
       end
       
       field_names_new{iObj} = sprintf('%s%s%s%s',cur_field_name,pre_str,mid_str,post_str);
    end
    
    [mask,index_new__value_old] = ismember(field_names_new,fields_old);
    
    %   field_names_new(~mask)'
    
    [mask2,loc2] = ismember(fields_old,field_names_new);
    
    %   fields_old(~mask2)'


end

function h__generateOldStats()

%NOTE: This function hasn't yet been exposed ...

    %======================================================================
    %======================================================================
    %                     Old code for testing
    %======================================================================
    %======================================================================
    %
    %   At this point the files that I need have been generated and this
    %   section doesn't need to be run ...
    %
    
    
    
    if false
    
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
    
    end
    
    
    if false
    %APPROACH 2  ==========================================================
    %Index exceeds matrix dimensions.
    
    expt_hist_file_old_v2 = fullfile(hist_path,'expt_hist_file_old_v2.mat');
    ctl_hist_file_old_v2  = fullfile(hist_path,'ctl_hist_file_old_v2.mat');
    
    worm2histogram(expt_hist_file_old_v2,expt_files);
    worm2histogram(ctl_hist_file_old_v2,ctl_files);
    
    stats_info_path_v2 = fullfile(hist_path,'stats_info_v2.mat');
    
    worm2StatsInfo(stats_info_path_v2, expt_hist_file_old_v2,[],[],ctl_hist_file_old_v2)
    
    end
    
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

end