classdef specs < handle
    %
    %   Class:
    %   seg_worm.stats.specs
    
    properties
       feature_field
       feature_category
       resolution
       is_zero_bin
       is_signed
       name
       short_name
       units
    end
    
    methods (Static,Hidden)
        function objs = getObjectsHelper(csv_path,class_function_handle,prop_names,prop_types)
           %
           %    The inherited objects can give relatively simple
           %    instructions on how their properties should be interpreted
           %    from their csv specification file.
           %
           %
           %    seg_worm.stats.specs.getObjectsHelper
           %
           %    INPUTS
           %    ===========================================================
           %
           %    TODO: Cleanup and finish documentation
           %
            %It would be nice to do the reading and object construction in 
           %here but Matlab is awkward for dynamic object creation 
           
           %TODO: It would be nice to improve this function to do the
           %casting inside this function ...
           output = sl.io.readDelimitedFile(csv_path,',',...
               'remove_empty_lines',true,'remove_lines_with_no_content',true);
           
           %Subtract header row ...
           n_objs = size(output,1) - 1;
           objs(n_objs) = class_function_handle();
           
           n_fields = length(prop_names);
           for iField = 1:n_fields
              %NOTE: we skip a header row :/
              cur_field_values = output(2:end,iField);
              cur_field_name   = prop_names{iField};
              cur_field_type   = prop_types(iField);
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
    methods
        function value = getLongField(obj)
           value = obj.feature_field;
        end
    end
    
end

