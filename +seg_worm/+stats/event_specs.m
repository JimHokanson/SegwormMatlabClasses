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
    
end

