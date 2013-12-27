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
    %   TODO:
    %   ================================================================
    %   1) Enumerate steps
    %   2) Create method for aggreating different files ...
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
    %
    %   JAH: I'm thinking that I'll move some of the old stats stuff to the
    %   hist class, and that this class will become explicitly for
    %   comparing the experiment and control data
    %
    %
    %
    
  
    
    
    properties
        
       %TODO: Move to object that both hist and stats display
       %Definitions in: seg_worm.stats.hist
       name
       short_name
       units
       feature_category
       hist_type
       motion_type
       data_type 
        
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
        function objs = stats(hist_objs)
            %
            %   seg_worm.stats()
            %

            if nargin == 0
                return
            end
            objs = seg_worm.stats.initObject(hist_objs);
        end
    end
    
    methods (Static)
       %seg_worm.stats.initObject(hist_objs)
       stats_objs = initObject(hist_objs)
    end
    
end

