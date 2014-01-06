classdef manager < handle
    %
    %   Class:
    %   seg_worm.stats.manager
    
    properties
    end
    
    methods
        function obj = manager(exp_hist_man,ctl_hist_man)
           
           obj.initObject(exp_hist_man.hists,ctl_hist_man.hists) 
            
            
        end
    end
    
end

