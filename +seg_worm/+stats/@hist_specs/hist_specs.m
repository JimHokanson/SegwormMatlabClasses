classdef hist_specs < handle
    %
    %   Class:
    %   seg_worm.stats.hist_specs
    %   
    %
    %   What is this class? Is it used??????
    %
    
    
    properties
       field            %...
       sub_fields       %... TODO: I think this should be singular ...
       %   I don't think this is ever plural
       type             %m ...
       category         %
       is_time_series   %
       resolution       %
       is_zero_bin      %
       is_signed        %
       name             %
       short_name       %
       unit             %
    end
    
    methods
        function obj = hist_specs(info1,info2)
            
           %TODO: Fix this ...
           if length(info2) > 1
               info2 = info2(1);
           end
            
           %NOTE: Eventually this will be replaced with
           obj.field          = info1.field;
           obj.sub_fields     = info1.subFields;
           obj.type           = info1.type;
           obj.category       = info1.category;
           obj.is_time_series = info1.isTimeSeries;
           obj.resolution     = info2.resolution;
           obj.is_zero_bin    = info2.isZeroBin;
           obj.is_signed      = info2.isSigned;
           obj.name           = info2.name;
           obj.short_name     = info2.shortName;
           obj.unit           = info2.unit;
           
        end
    end
    
end

