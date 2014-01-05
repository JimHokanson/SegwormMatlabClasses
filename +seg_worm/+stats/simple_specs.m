classdef simple_specs < seg_worm.stats.specs
    %
    %   Class:
    %   seg_worm.stats.simple_specs
    %
    
    properties
    end
    
    methods (Static)
        function objs = getSpecs()
            %
            %
            %   s_specs = seg_worm.stats.simple_specs.getSpecs();
            %
            %
            
           csv_path = fullfile(sl.dir.getMyBasePath,'docs','simple_features.csv');
           
           info = {'feature_field'      1
                   'feature_category'   1
                   'resolution'         2
                   'name'               1
                   'short_name'         1
                   'units'              1};
           
           
           prop_names = info(:,1);
           prop_types = cat(1,info{:,2});

           fh = @seg_worm.stats.simple_specs;
           
           objs = seg_worm.stats.specs.getObjectsHelper(csv_path,fh,prop_names,prop_types);
            
           [objs(:).is_zero_bin] = deal(false);
           [objs(:).is_signed]   = deal(false);
        end
    end
    methods
        function data = getData(obj,feature_obj)
           data  = eval(['feature_obj.' obj.feature_field]); 
        end 
    end
    
end

