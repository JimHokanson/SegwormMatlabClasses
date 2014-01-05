classdef stats < handle
    %
    %   Class:
    %   seg_worm.stats
    %
    %   Important Files
    %   ---------------------------------------------------------
    %   seg_worm.w.stats.worm2histogram         TODO: Move this ...
    %   seg_worm.w.stats.addWormHistograms
    %   seg_worm.w.stats.worm2StatsInfo
    %   seg_worm.w.stats.wormStats2Matrix
    %
    %
    %   QUESTIONS:
    %   ================================================================
    %   1)
    %
    %   Procesing Steps:
    %   ================================================================
    %   1) Create histogram files for "experiment" and "control" files
    %       - seg_worm.stats.hist
    %   2)
    %
    %
    %   Some of the statistics are aggegrate:
    %   p_value
    %   q_value
    %   list of exclusive features
    
  
    
    
    properties
        
       %TODO: Move to object that both hist and stats display
       %
       %ALSO: We need two, one for experiment and one for controls
       %Definitions in: seg_worm.stats.hist
       name
       short_name
       units
       feature_category
       hist_type
       motion_type
       data_type 
        
       %New properties
       %-------------------------------------------------------------------
       q_normal_experiment
       q_normal_control
       z_score_experiment
       z_score_control
       p_significance
       p_normal
       
       %-------------------------------------------------------------------
       z_score   %not populated if no controls are provided ...
       mean      %mean of the mean hist values
       std       %std of the hist values
       n_samples %# of videos where the mean is not NaN
       p_normal = NaN  %probability of being a normal distribution
       %
       %    seg_worm.fex.swtest(data(i).dataMeans, 0.05, 0)
       q_normal  %
    end
    
    methods
        function obj = stats()
            %
            %   seg_worm.stats()
            %
            
            %obj = obj.initObject(exp_hist,ctl_hist);
        end
    end
    
    methods (Static)
       %seg_worm.stats.initObject(hist_objs)
       stats_objs = initObject(hist_objs)
    end
    
end

