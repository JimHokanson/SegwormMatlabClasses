function data = catWormData(wormFile, varargin)
%catWormData   Concatenate the worm data.
%
%   data = seg_worm.w.util.catWormData(wormFile, index, start_frame, end_frame)
%   
%
%   NOT CALLED BY ANY FUNCTIONS IN CODE BASE
%
%
%
%   DATA = CATWORMDATA(WORMFILE)
%   DATA = CATWORMDATA(WORMFILE, INDEX)
%   DATA = CATWORMDATA(WORMFILE, INDEX, STARTFRAME)
%   DATA = CATWORMDATA(WORMFILE, INDEX, STARTFRAME, ENDFRAME)
%
%   Inputs:
%       See WORM_FILE in OldWormFormats.m
%
%       index      - the data index in the block;
%                    if empty, return all data
%       startFrame - the first frame to use;
%                    if empty, we start at the first frame
%       endFrame   - the last frame to use;
%                    if empty, we end at the last frame
%
%   Output:
%
%       data - the concatenated worm data
%
%   See also:
%   SAVEWORMVIDEOFRAMES
%   NORMWORMS
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Determine the data index.
index = 1:12;
if ~isempty(varargin)
    index = varargin{1};
end

% Determine the start and end frames.
startFrame = [];
if length(varargin) > 1
    startFrame = varargin{2};
end
endFrame = [];
if length(varargin) > 2
    endFrame = varargin{3};
end

% Get the data blocks.
state = index;
dataBlocks = worm2func(@echoFunc, state, wormFile, startFrame, ...
    endFrame, 0, 0);

% Concatenate the data blocks.
data = cell(length(dataBlocks{1}), 1);
for i = 1:length(dataBlocks)
    for j = 1:length(data)
        if ~isempty(dataBlocks{i})
            switch ndims(dataBlocks{i}{j})
                case 2
                    data{j} = cat(2, data{j}, dataBlocks{i}{j});
                case 3
                    data{j} = cat(3, data{j}, dataBlocks{i}{j});
                otherwise
                    error('catWormBlocks:BadVariable', ['Data cell ' ...
                        num2str(i) ' has inappropropriate dimensionality']);
            end
        end
    end
end
end

% Echo the data.
function [data,dataBlockI] = echoFunc(dataInfo, dataBlockI)

% Initialize the variables.
dataBlock = dataInfo.data;
startI    = dataInfo.startDataFrameI;
endI      = dataInfo.endDataFrameI;

% Extract the subset of the data.
data = cell(length(dataBlockI), 1);
for i = 1:length(data)
    j = dataBlockI(i);
    switch ndims(dataBlock{j})
        case 2
            data{i} = dataBlock{j}(:,startI:endI);
        case 3
            data{i} = dataBlock{j}(:,:,startI:endI);
        otherwise
            error('catWormBlocks:BadVariable', ['Data cell ' ...
                num2str(j) ' has inappropropriate dimensionality']);
    end
end
end
