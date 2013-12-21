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
       feature_field
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
           
           %TODO: It would be nice to improve this function to do the
           %casting inside this function ...
           output = sl.io.readDelimitedFile(csv_path,',');
           
           %TODO: These properties would be better if paired then split ...
           %i.e. group names and type, then split
           %
           %    info = {'feature_field' 1
           %            'index' 2 
           %            'feature_category' 1, etc
           %
           %These are the property names that we will assign each column to
           field_names = {'feature_field' 'index' 'feature_category' 'is_time_series' ...
               'resolution' 'is_zero_bin' 'is_signed' 'name' 'short_name' 'units'};
           
           %1 - strings
           %2 - numeric
           %3 - logical
           field_type = [1 2 1 3 2 3 3 1 1 1];
           
           %Subtract header row ...
           n_objs = size(output,1) - 1;
           
           objs(n_objs) = seg_worm.stats.movement_specs;
           
           n_fields = length(field_names);
           for iField = 1:n_fields
              %NOTE: we skip a header row :/
              cur_field_values = output(2:end,iField);
              cur_field_name   = field_names{iField};
              cur_field_type   = field_type(iField);
              switch cur_field_type
                  case 1
                      %strings, do nothing ...
                  case 2
                      cur_field_values = num2cell(str2double(cur_field_values));
                  case 3
                      cur_field_values = num2cell(cellfun(@(x) x == '1',cur_field_values));
              end
              [objs.(cur_field_name)] = deal(cur_field_values{:});
                      
           end
          
        end
    end
    
end

