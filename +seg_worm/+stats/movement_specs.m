classdef movement_specs < seg_worm.stats.specs
    %
    %   Class:
    %   seg_worm.stats.movement_specs
    %
    %   This class specifies how to treat each movement related feature for
    %   histogram processing.
    %
    %
    %   Access via static method:
    %   seg_worm.stats.movement_specs.getSpecs()
    %
    %   See Also:
    %   seg_worm.stats.hist.createHistograms
    %
    %   TODO:
    %   - might need to incorporate seg_worm.w.stats.wormStatsInfo
    %   - remove is_time_series entry ...
    
    properties
       index
       %feature_category
       is_time_series %TODO: This can be removed ...
       %resolution
       %is_zero_bin %This might not be important
       %is_signed   %I think this dictates having 4 or 16 events ...
%        name
%        short_name
%        units
    end

    methods (Static)
        function objs = getSpecs()
           %seg_worm.stats.movement_specs.getSpecs();
           
           csv_path = fullfile(sl.dir.getMyBasePath,'movement_features.csv');
           

           
           %TODO: These properties would be better if paired then split ...
           %i.e. group names and type, then split
           %
           %    info = {'feature_field' 1
           %            'index' 2 
           %            'feature_category' 1, etc
           %
           %These are the property names that we will assign each column to
           prop_names = {'feature_field' 'index' 'feature_category' 'is_time_series' ...
               'resolution' 'is_zero_bin' 'is_signed' 'name' 'short_name' 'units'};
           
           %1 - strings
           %2 - numeric
           %3 - logical
           prop_types = [1 2 1 3 2 3 3 1 1 1];
           
           fh = @seg_worm.stats.movement_specs;
           
           objs = seg_worm.stats.specs.getObjectsHelper(csv_path,fh,prop_names,prop_types);

           

          
        end
    end
    
end

