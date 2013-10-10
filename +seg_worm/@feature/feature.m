classdef feature
    %
    %   Class:
    %   seg_worm.feature
    %
    %   See Also:
    %   seg_worm.feature.roots
    %   seg_worm.feature.displayInfo
    
    
    %????
    %   - How should I delinieate between the features in the file
    %   and the features that are used for statistics ...
    
    
    %??? What really breaks up the different feature types???
    %
    %   - the type?
    
%From seg_worm.feature.displayInfo
%------------------------------------------------------------------
%              name       = a descriptive name
%              shortName  = a shortened name
%
%              And each feature has the following subfields:
%
%              resolution = the resolution for the histogram
%              isZeroBin  = is the histogram centered at zero?
%              isSigned   = is the data signed (+/-)?
%              name       = a descriptive feature name (for titles)
%              shortName  = a shortened feature name (for legends)
%              unit       = the feature's measurement unit (for labels)
    
    
    
    
    properties
       field
       sub_fields
       type
       category
       is_time_series
    end
    
    methods
    end
    
end

