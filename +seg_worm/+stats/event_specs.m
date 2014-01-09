classdef event_specs < seg_worm.stats.specs
    %
    %   Class:
    %   seg_worm.stats.event_specs
    %
    
    properties
        sub_field
    end
    
    methods (Static)
        function objs = getSpecs()
            %
            %
            %   e_specs = seg_worm.stats.event_specs.getSpecs();
            %
            %
            
            csv_path = fullfile(sl.dir.getMyBasePath,'docs','event_features.csv');
            
            info = {'feature_field'      1
                'sub_field'          1
                'feature_category'   1
                'resolution'         2
                'is_signed'          3
                'name'               1
                'short_name'         1
                'units'              1};
            
            
            prop_names = info(:,1);
            prop_types = cat(1,info{:,2});
            
            fh = @seg_worm.stats.event_specs;
            
            objs = seg_worm.stats.specs.getObjectsHelper(csv_path,fh,prop_names,prop_types);
            
            [objs(:).is_zero_bin] = deal(false);
        end
    end
    methods
        function data = getData(obj,feature_obj)
            %NOTE: Because we are doing structure array indexing, we need to capture
            %multiple outputs using [], otherwise we will only get the first value
            %...
            data  = eval(['feature_obj.' obj.feature_field]);
            
            if ~isempty(data) && ~isempty(obj.sub_field)
                %This will go from:
                %   frames (structure array)
                %to:
                %   frames.time
                %for example.
                %
                %It is also used for event.ratio.time and event.ratio.distance
                %     going from:
                %         ratio (structure or [])
                %     to:
                %         ratio.time
                %         ratio.distance
                %
                %
                data = [data.(obj.sub_field)];
            end
        end
        function value = getLongField(obj)
           value = obj.feature_field;
           if ~isempty(obj.sub_field)
              value = [value '.' obj.sub_field]; 
           end
        end
    end
    
end

