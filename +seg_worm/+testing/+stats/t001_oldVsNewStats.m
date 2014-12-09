function t001_oldVsNewStats
    % To run the old code the folder "oldStats" needs to be on the path:
    addpath('C:\Users\mcurrie\Desktop\GitHub\SegwormMatlabClasses\oldStats')

    %NOTE: These files don't necessarily make a valid set of comparisons,
    %we are just grabbing a set of feature files, not a set that reproduces
    %a specific worm.
    %root_path = 'C:\Users\RNEL\Google Drive\open_worm\example_data\30m_wait';
    root_path = 'C:\Backup\Google Drive\OpenWorm\OpenWorm Public\movement_validation\example_data\30m_wait';

    %This should be an empty folder for saving things into temporarlly for
    %multiple runs since some of the old code is quite slow
    %temp_save_path = 'F:\worm_data\stats_testing';
    temp_save_path = 'C:\Backup\stats_testing';

    old_vs_new_stats(root_path, temp_save_path)
end


function old_vs_new_stats(root_path, temp_save_path)
%
%   Compares old vs. new code for statistics.  As a result, we need
%   to use the original code from MRC in addition to Jim's modified 
%   versions which have been cleaned up (the "new" code).
%
%   Feature Count
%   --------------------------------------------------
%   726 - old files
%   708 - new files (events need to be signed)
%   702 - nature methods (removal of paused body motion, 6 categories, 4
%   features each paused all, paused abs., paused neg, paused pos, 24
%   total)
%
%   Current Status:
%   ------------------------------------------------------------------
%   1) This code is really meant to be run via evaluation. It needs a bit
%   of cleaning up to run as a function. It also needs to have test results
%   actually by displayed
%   2) It would be nice if the test files were passed in.
%   3) The new method and old method seem to match.
%

    feature_files = sl.dir.rdir([root_path '\**\*.mat']);
    
    %We'll take the first 10 for the "experiment" and the next 10 for the
    %"control"
    %expt_files = {feature_files(1:3).name};
    %ctl_files  = {feature_files(14:16).name};
    expt_files = {feature_files(1:10).name};
    %ctl_files  = {feature_files(14:23).name};
    ctl_files  = {feature_files(11:20).name};

    
    
    % STEP 1: Compute Old Histograms - write to disk
    %----------------------------------------------------------------------
    sl.dir.createFolderIfNoExist(temp_save_path);
    [stats_info_path_v3, stats_matrix_path_v3] ...
        = compute_old_histograms(temp_save_path, expt_files, ctl_files);

    
    % STEP 2: Compute "new" Histograms - (just save in memory, don't 
    %                                    bother writing to disk)
    %----------------------------------------------------------------------    
    % TODO: Make this seg_worm.getHists for short
    hist_man_exp = seg_worm.stats.hist.manager(expt_files);
    hist_man_ctl = seg_worm.stats.hist.manager(ctl_files);
    % TODO: Make this seg_worm.getStats for short
    stats_manager = seg_worm.stats.manager(hist_man_exp,hist_man_ctl);

    
    % STEP 3: Compare old vs. new statistics, display the results
    %----------------------------------------------------------------------
    compare_old_vs_new_stats(stats_info_path_v3, ...
                             stats_matrix_path_v3, ...
                             hist_man_exp, hist_man_ctl, stats_manager)
    
end



function compare_old_vs_new_stats(stats_info_path_v3, ...
                                  stats_matrix_path_v3, ...
                                  hist_man_exp, hist_man_ctl,stats_manager)
%==========================================================================    
%                   Comparison Code    
%========================================================================== 
    
    h_si    = load(stats_info_path_v3);   %stats info   - old objects
    h_sm    = load(stats_matrix_path_v3); %stats matrix - old objects
    h_stats = stats_manager.stats;        %               new objects

    [index_new__value_old, field_names_new, fields_old] = ...
        h_getOrder(h_si,h_sm,h_stats);


    %h_si - contains    


    field_names_new_ordered(index_new__value_old) = field_names_new;
    h_new(index_new__value_old) = h_stats;

    clear h_exp_new
    clear h_ctl_new
    h_exp_new(index_new__value_old) = hist_man_exp.hists;
    h_ctl_new(index_new__value_old) = hist_man_ctl.hists;

    %??? - why were some of the names empty????
    keep_mask    = ~cellfun('isempty',{h_new.name});
    kept_indices = find(keep_mask);

    %We don't have the mean or std data here
    %mean control
    %mean_control_old = h_sm.control.stats.mean(keep_mask);
    %mean_control_new = [



    %zscore comparision
    %--------------------------------------------------------------------------
    %zscore_control - 0 by definition
    %
    %zscore_experiment - looks the same except for NaN mismatch due to 
    %differences in definition ...
    %--------------------------------------------------------------------------
    %old code - set in worm2StatsInfo
    fprintf('Zscore comparison --------------------------\n')
    zscore_experiment_old = h_sm.worm.stats.zScore(keep_mask);
    zscore_experiment_new = [h_new(keep_mask).z_score_experiment];

    %difference in the z_scores, 1 is a bit high with the 30m_wait set, but
    %isn't too bad (0.01)
    d_zscore_experiment = zscore_experiment_old - zscore_experiment_new;
    %find fields where only 1 is nan
    nan_mismatch_zscore = find(xor(isnan(zscore_experiment_old),isnan(zscore_experiment_new)));

    %These are being displayed but they probably don't need to be
    %
    %Many olds are NaN, not -Inf or Inf as the documentation states 
    zscore_experiment_old(nan_mismatch_zscore)
    zscore_experiment_new(nan_mismatch_zscore)
    fields_old(nan_mismatch_zscore)'

    %TODO: What to summarize here for z_score




    %p_t comparison - values same, except for exclusive values ...
    %--------------------------------------------------------------------------
    p_t_old = h_sm.worm.sig.pTValue;
    p_t_new = [h_new.p_t];

    d_p_t = p_t_old - p_t_new;

    nan_mismatch_p_t = find(xor(isnan(p_t_old),isnan(p_t_new)));
    %433   571   572   573   574   587   588   589   590


    %p_w comparison - values look the same
    %-------------------------------------------------------------------
    p_w_old = h_sm.worm.sig.pWValue;
    p_w_new = [h_new.p_w];

    d_p_w = p_w_old - p_w_new;

    nan_mismatch_p_w = find(xor(isnan(p_w_old),isnan(p_w_new)));


    %I think the differences might be because of the NaN mismatch between
    %between the p-values, since the q values are dependent on the aggregate.

    %q_t comparison - values look the same ...
    %-------------------------------------------------------
    q_t_old = [h_si.significance.features.qTValue];
    %q_t_old = h_sm.worm.sig.qTValue.all;
    %q_t_old = h_sm.worm.sig.qTValue.strain;
    q_t_new = [h_new.q_t];

    d_q_t = q_t_old - q_t_new;

    nan_mismatch_q_t = find(xor(isnan(q_t_old),isnan(q_t_new)));

    %q_w comparison - strain is relatively close ...
    %-------------------------------------------------------------
    q_w_old = h_sm.worm.sig.qWValue.all;
    q_w_old = h_sm.worm.sig.qWValue.strain;
    q_w_new = [h_new.q_w];

    d_q_w = q_w_old - q_w_new;

    nan_mismatch_q_w = find(xor(isnan(q_w_old),isnan(q_w_new)));

end

function [stats_info_path_v3, stats_matrix_path_v3] = ...
        compute_old_histograms(temp_save_path, expt_files, ctl_files)
    %Step 1: Compute Old Histograms - write to disk
    %----------------------------------------------------------------------
    old_hist_save_path =  sl.dir.createFolderIfNoExist(temp_save_path,'histograms');
    
    %This process is slow - many minutes, we'll save the results to disk
    %for faster running later. Additionally, the old files are built to
    %handle reading and writing things from disk, not from memory
    
    all_expt_hists_v3 = cell(1,10);
    for iFile = 1:10
        cur_file_path   = expt_files{iFile};
        [~,cur_name]    = fileparts(cur_file_path);
        cur_output_path = fullfile(old_hist_save_path,[cur_name '_v3.mat']);
        %If this is not on the path, you need to add the folder 'oldStats'
        if ~exist(cur_output_path,'file')
            worm2histogram(cur_output_path, cur_file_path)
        end
        all_expt_hists_v3{iFile} = cur_output_path;
    end
    
    all_ctl_hists_v3 = cell(1,10);
    for iFile = 1:10
        cur_file_path   = ctl_files{iFile};
        [~,cur_name]    = fileparts(cur_file_path);
        cur_output_path = fullfile(old_hist_save_path,[cur_name '_v3.mat']);
        if ~exist(cur_output_path,'file')
            worm2histogram(cur_output_path, cur_file_path)
        end
        all_ctl_hists_v3{iFile} = cur_output_path;
    end

    %We've previously computed individual histogram files for each video.
    %The stats are computed on a set of hisograms versus another set. To
    %compute the merged set, the following code is run.
    expt_hist_merged = fullfile(old_hist_save_path,'expt_hist_file_old_merged_v3.mat');
    ctl_hist_merged  = fullfile(old_hist_save_path,'ctl_hist_file_old_merged_v3.mat');
    
    if ~exist(expt_hist_merged,'file')
        addWormHistograms(expt_hist_merged,all_expt_hists_v3);
    end
    if ~exist(ctl_hist_merged,'file')
        addWormHistograms(ctl_hist_merged,all_ctl_hists_v3);
    end
    
    %Compute old stats
    % Take the 10 experiments and 10 controls and compares them
    % e.g. mean fwd velocity per video across ten videos
    %      so take those 10 and compare to the other 10 and 
    %      calculates p's and q's from the null hypothesis that 
    %      the figures were drawn from the same distribution
    %-----------------
    stats_info_path_v3 = fullfile(temp_save_path,'stats_info_v3.mat');
    stats_matrix_path_v3 = fullfile(temp_save_path,'stats_matrix_v3.mat');
    if ~exist(stats_matrix_path_v3,'file')
        worm2StatsInfo(stats_info_path_v3, expt_hist_merged,[],[],ctl_hist_merged)
        wormStats2Matrix(stats_matrix_path_v3, stats_info_path_v3)
    end

end

function [index_new__value_old,field_names_new,fields_old] = h_getOrder(h_si,h_sm,h_stats)
%
%
%   This function is a hack to align the old and new stats ordering (look
%   away :/ )
%
%


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


    n_old = length(h_si.controlData);

    fields_old     = cell(1,n_old);
    categories_all = cell(1,n_old);
    type_all       = cell(1,n_old);
    sub_type_all   = cell(1,n_old);
    sign_all       = cell(1,n_old);
    for iObj = 1:n_old
       fields_old{iObj}     = h_si.dataInfo(iObj).field.histogram;
       categories_all{iObj} = h_si.dataInfo(iObj).category;
       type_all{iObj}       = h_si.dataInfo(iObj).type;
       sub_type_all{iObj}   = h_si.dataInfo(iObj).subType;
       sign_all{iObj}       = h_si.dataInfo(iObj).sign;
    end
        
    n_new = length(h_stats);
    
    field_names_new = cell(1,n_new);
    for iObj = 1:n_new       
       cur_stat = h_stats(iObj);
       cur_field_name = cur_stat.field;

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
       
       if strcmp(cur_stat.hist_type,'event')
           
           %Different rules for frames and no frames
           
           I = strfind(cur_field_name,'frames');
           
           if ~isempty(I)
           
               cur_field_name = regexprep(cur_field_name,'\.frames\.','\.');
               pre_str  = '';
               mid_str = '.histogram.data.mean';
           else
                pre_str  = '';
                mid_str = '.data';
                post_str = '';
           end
       elseif strcmp(cur_stat.hist_type,'simple')
           cur_field_name = regexprep(cur_field_name,'\.times','');
           pre_str  = '';
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

           mid_str = '.histogram.data.mean';
       end
       
       field_names_new{iObj} = sprintf('%s%s%s%s',cur_field_name,pre_str,mid_str,post_str);
    end
    
    
    %NOTE: This needs the better ismember function that I've been meaning
    %to right for a while ...
    
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