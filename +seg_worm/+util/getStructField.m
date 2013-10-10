function data = getStructField(data, field, varargin)
%getStructField  Get a field in a structure.
%
%   data = getStructField(DATA, FIELD)
%
%   data = getStructField(data, field, *isSafe)
%
%   seg_worm.util.getStructField
%
%
%
%   NOTE: If a field does not exist then data will be returned as empty (no
%   error)
%
%
%   Inputs:
%       data   - the structure
%       field  - the field to get (a period delimited string)
%       isSafe - are we using safe mode? when using safe mode, each field
%                is tested to ensure it exists; the default is no (false)
%
%
%
%
%   Output:
%       data - the field data;
%              empty if the field cannot be found
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.



% Are we using safe mode?
isSafe = false;
if ~isempty(varargin)
    isSafe = varargin{1};
end


if isempty(field)
    % No data.
    
    %Why would field 'input' be empty ?????
    %I could understand the field not existing, but the request being
    %empty doesn't make sense to me
    data = [];
    return;
else
    % Get the field.
    fields = textscan(field, '%s', 'Delimiter', '.');
    fields = fields{1};
    
    % Check each field to ensure it exists before traversing.
    if isSafe
        for i = 1:length(fields)
            
            % Have we reached the end?
            if isempty(data) || ~isstruct(data) || ~isfield(data, fields{i})
                data = [];
                return;
            end
            
            % Go deeper.
            data = [data.(fields{i})];
        end
        
    % Traverse each field without checking to ensure it exists.
    else
        for i = 1:length(fields)
            
            % Have we reached the end?
            if isempty(data) || ~isstruct(data)
                data = [];
                return;
            end
            
            % Go deeper.
            data = [data.(fields{i})];
        end
    end
end
end
