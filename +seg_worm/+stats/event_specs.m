classdef event_specs < seg_worm.stats.specs
    %
    %   Class:
    %   seg_worm.stats.event_specs
    %
    
    properties
        sub_field
        signed_field = '' %True will indicate that the data should be negated ...
        make_zero_if_empty
        remove_partials
    end
    
    methods (Static)
        function objs = getSpecs()
            %
            %
            %   e_specs = seg_worm.stats.event_specs.getSpecs();
            %
            %
            
            csv_path = fullfile(sl.stack.getMyBasePath,'docs','event_features.csv');
            
            %See seg_worm.stats.specs.getObjectsHelper for details
            info = {...
                'feature_field'         1
                'sub_field'             1
                'feature_category'      1
                'resolution'            2
                'signed_field'          1
                'name'                  1
                'short_name'            1
                'units'                 1
                'make_zero_if_empty'    3
                'remove_partials'       3};
            
            
            prop_names = info(:,1);
            prop_types = cat(1,info{:,2});
            
            fh   = @seg_worm.stats.event_specs;
            
            objs = seg_worm.stats.specs.getObjectsHelper(csv_path,fh,prop_names,prop_types);
            
            %Finish populating inherited methods
            %-------------------------------------------------
            [objs(:).is_zero_bin] = deal(false);
            
            signed_field_names  = {objs.signed_field};
            has_name            = ~cellfun('isempty',signed_field_names);
            [objs(:).is_signed] = sl.struct.dealArray(has_name);
        end
    end
    methods
        function data = getData(obj,feature_obj,n_samples)
            %NOTE: Because we are doing structure array indexing, we need to capture
            %multiple outputs using [], otherwise we will only get the first value
            %...
            
            %TODO: Perhaps we'll add a property in the new object instead
            %of this poor check ...
            is_old_code = isstruct(feature_obj);

            if is_old_code
                start_value = 0;
                end_value   = n_samples; %BUG IN OLD CODE: n_samples matches
                %behavior, n_samples -1 does not
            else
                start_value = 1;
                end_value   = n_samples;
            end

            data = sl.struct.getSubField(feature_obj,obj.feature_field);
            
            if ~isempty(data)
                
                if ~isempty(obj.sub_field)
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
                    parent_data = data;
                    
                    data = [data.(obj.sub_field)];
                    
                    if obj.is_signed
                        negate_mask = [parent_data.(obj.signed_field)];
                        data(negate_mask) = -1*data(negate_mask);
                    end
                    
                    if obj.remove_partials
                       starts = [parent_data.start];
                        ends   = [parent_data.end];

                        remove_mask = false(1,length(starts));

                        if starts(1) == start_value
                            remove_mask(1) = true;
                        end

                        if ends(end) == end_value
                            remove_mask(end) = true;
                        end 
                        
                        data(remove_mask) = [];
                    end
                    
                else
                    %Check things that don't currently make sense unless
                    %nested in the way that we expect (i.e. in the frames
                    %struct)
                    
                    %TODO: Can't be signed
                    %TODO: Can't remove partials
                    
                end
                
            end
            
            if isempty(data) && obj.make_zero_if_empty
                data = 0;
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

