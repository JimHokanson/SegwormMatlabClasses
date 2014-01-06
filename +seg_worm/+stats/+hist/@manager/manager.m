classdef manager < handle
    %
    %   Class:
    %   seg_worm.stats.hist.manager
    
    properties
       %This should eventually hold the video information ...
       %.info [n_videos x 1]
       
       hists
    end
    
    methods
        function obj = manager(feature_objs_or_paths)
        
           
        %Loop over all feature files and get histogram objects for each
        %--------------------------------------------------------------------------
        if ischar(feature_objs_or_paths)
            feature_objs_or_paths = {feature_objs_or_paths};
        end

        n_videos = length(feature_objs_or_paths);
        hist_cell_array = cell(n_videos,1);

        for iVideo  = 1:n_videos
            cur_path_or_obj = feature_objs_or_paths{iVideo};
            if ischar(cur_path_or_obj)
                h = load(cur_path_or_obj);
                feature_obj = h.worm;
            else
                feature_obj = cur_path_or_obj;
            end
            
            %TODO: Need to add on info to properties 
            %feature_obj.info -> obj.info
            
            hist_cell_array{iVideo} = initObjects(obj,feature_obj);
        end

        %Merge the objects from each file
        %--------------------------------------------------------------------------
        obj.hists = seg_worm.stats.hist.mergeObjects(hist_cell_array); 
            
        end
        %TODO: Should also have a merge manager objects function as well
        %...
    end
    
end

