function bends = getLocomotionBends(wormFile, varargin)
%WORMBENDS Compute the temporal bending frequency at the nose, head,
%midbody, and tail.
%
%
%   bends =
%   seg_worm.feature_helpers.locomotion.getLocomotionBends(worm_file,motion_mode,ventral_mode)
%
%
%   Old Name: wormBends.m
%
%   Unfinished  Unfinished  Unfinished  Unfinished
%
%       motionMode  - the locomotion mode. An optional argument, that when
%                     present, removes non-foraging bending frequencies
%                     measured during pauses in locomotion (e.g.,
%                     frequencies measured as a result of defecation).
%                     The modes are:
%
%                     -1 = backward locomotion
%                      0 = no locomotion (the worm is paused)
%                      1 = forward locomotion
%
%       ventralMode - the ventral side mode:
%
%                     0 = the ventral side is unknown
%                     1 = the ventral side is clockwise
%                     2 = the ventral side is anticlockwise
%
%   Outputs:
%       bends - the bend information as a struct with the subfields for the
%               "nose", "head", "midbody", and "tail" bends; each of these 
%               subfields is a struct containing their "amplitude" and
%               "frequency", except foraging which contains an "amplitude"
%               and an "angleSpeed" (an angular speed in place of the
%               frequency). Each value maintains its dorsal/ventral sign.
%
%
%   Nature Methods Description - Crawling and Foraging (might split)
%   =======================================================================
%
%   Crawling. 
%
%   Worm crawling is expressed as both an amplitude and frequency
%   (Supplementary Fig. 4e). We measure these features instantaneously at
%   the head, midbody, and tail. The amplitude and frequency are signed
%   negatively whenever the worm’s ventral side is contained within the
%   concave portion of its instantaneous bend.
%  
%   Crawling is only measured during forward and backward motion states.
%   The worm bend mean angles (described in the section on “Posture”) show
%   a roughly periodic signal as the crawling wave travels along the worm’s
%   body. This wave can be asymmetric due to differences in dorsal-ventral
%   flexibility or simply because the worm is executing a turn. Moreover
%   the wave dynamics can change abruptly to speed up or slow down.
%   Therefore, the signal is only roughly periodic and we measure its
%   instantaneous properties.
%
%   Worm bends are linearly interpolated across unsegmented frames. The
%   motion states criteria (described earlier in this section) guarantee
%   that interpolation is no more than 1/4 of a second long. For each
%   frame, we search both backwards and forwards for a zero crossing in the
%   bend angle mean – the location where the measured body part (head,
%   midbody, or tail) must have hit a flat posture (a supplementary bend
%   angle of 0°). This guarantees that we are observing half a cycle for
%   the waveform. Crawling is bounded between 1/30Hz (a very slow wave that
%   would not resemble crawling) and 1Hz (an impossibly fast wave on agar).
%   If the window between zero crossings is too small, the nearest zero
%   crossing is assumed to be noise and we search for the next available
%   zero crossing in its respective direction. If the window is too big,
%   crawling is marked undefined at the frame. Once an appropriate window
%   has been found, the window is extended in order to center the frame and
%   measure instantaneous crawling by ensuring that the distance on either
%   side to respective zero crossings is identical. If the distances are
%   not identical, the distance of the larger side is used in place of the
%   zero-crossing distance of the smaller side in order to expand the small
%   side and achieve a symmetric window, centered at the frame of interest.
% 
%   We use a Fourier transform to measure the amplitude and frequency
%   within the window described above. The largest peak within the
%   transform is chosen for the crawling amplitude and frequency. If the
%   troughs on either side of the peak exceed 1/2 its height, the peak is
%   rejected for being unclear and crawling is marked as undefined at the
%   frame. Similarly, if the integral between the troughs is less than half
%   the total integral, the peak is rejected for being weak.
%
%   Foraging
%   
%   Foraging. Worm foraging is expressed as both an amplitude and an
%   angular speed (Supplementary Fig. 4g). Foraging is signed negatively
%   whenever it is oriented towards the ventral side. In other words, if
%   the nose is bent ventrally, the amplitude is signed negatively.
%   Similarly, if the nose is moving ventrally, the angular speed is signed
%   negatively. As a result, the amplitude and angular speed share the same
%   sign roughly only half the time. Foraging is an ambiguous term in
%   previous literature, encompassing both fine movements of the nose as
%   well as larger swings associated with the head. Empirically we have
%   observed that the nose movements are aperiodic while the head swings
%   have periodicity. Therefore, we measure the aperiodic nose movements
%   and term these foraging whereas the head swings are referred to as
%   measures of head crawling (described earlier in this section).
% 
%   Foraging movements can exceed 6Hz7 and, at 20-30fps, our video frame
%   rates are just high enough to resolve the fastest movements. By
%   contrast, the slowest foraging movements are simply a continuation of
%   the crawling wave and present similar bounds on their dynamics.
%   Therefore, we bound foraging between 1/30Hz (the lower bound used for
%   crawling) and 10Hz.
% 
%   To measure foraging, we split the head in two (skeleton points 1-4 and
%   5-8) and measure the angle between these sections. To do so, we measure
%   the mean of the angle between subsequent skeleton points along each
%   section, in the tail-to-head direction. The foraging angle is the
%   difference between the mean of the angles of both sections. In other
%   words, the foraging angle is simply the bend at the head.
% 
%   Missing frames are linearly interpolated, per each skeleton point, for
%   fragments up to 0.2 seconds long (4-6 frames at 20-30fps – twice the
%   upper foraging bound). When larger fragments are missing, foraging is
%   marked undefined. Segmentation of the head at very small time scales
%   can be noisy. Therefore, we smooth the foraging angles by convolving
%   with a Gaussian filter 1/5 of a second long (for similar reasons to
%   those mentioned in frame interpolation), with a width defined by the
%   Matlab “gausswin” function’s default ? of 2.5 and normalized such that
%   the filter integrates to 1.
% 
%   The foraging amplitude is defined as the largest foraging angle
%   measured, prior to crossing 0°. In other words, the largest nose bend
%   prior to returning to a straight, unbent position. Therefore, the
%   foraging amplitude time series follows a discrete, stair-step pattern.
%   The amplitude is signed negatively whenever the nose points towards the
%   worm’s ventral side. The foraging angular speed is measured as the
%   foraging angle difference between subsequent frames divided by the time
%   between these frames. To center the foraging angular speed at the frame
%   of interest and eliminate noise, each frame is assigned the mean of the
%   angular speed computed between the previous frame and itself and
%   between itself and the next frame. The angular speed is signed
%   negatively whenever its vector points towards the worm’s ventral side.

% Are we using the locomotion mode to remove non-foraging bends?
pausedMotion = [];
if ~isempty(varargin)
    motionMode = varargin{1};
    pausedMotion = (motionMode == 0);
end

% Where is the ventral side located?
ventralMode = 0;
if length(varargin) > 1
    ventralMode = varargin{2};
end

% Check the worm file.
if ~exist(wormFile, 'file')
    error('multiScaleWorm:BadWormFile', ['Cannot find ''' wormFile '''']);
end

% Are the normalized blocks separate?
global blockFilePath;
blockFilePath = [];
load(wormFile, 'block1');
if ~exist('block1', 'var')
    
    % Get the path from the file.
    blockFilePath = fileparts(wormFile);
    
    % Use the current path.
    if isempty(blockFilePath)
        blockFilePath = pwd;
    end
    
% Clean up memory.
else
    clear('block1');
end

% Determine the variables (the blocks are in separate files).
if ~isempty(blockFilePath)
    
    % Load the worm information.
    load(wormFile, 'myAviInfo', 'normBlockList');
    
    % Check the variables.
    varName = 'myAviInfo';
    if ~exist(varName, 'var')
        error('worm2func:BadVariable', ...
            ['Cannot load ''' varName ''' from ''' wormFile '''']);
    else
        fps = myAviInfo.fps;
    end
    varName = 'normBlockList';
    if ~exist(varName, 'var')
        error('worm2func:BadVariable', ...
            ['Cannot load ''' varName ''' from ''' wormFile '''']);
    end
    
    % Determine the number of blocks.
    blocks = length(normBlockList);
    
    % Determine the block size.
    load(fullfile(blockFilePath, normBlockList{1}));
    eval(['block = ' normBlockList{1} ';']);
    clear(normBlockList{1});
    blockSize = size(block{1}, 2);
    
    % Determine the last frame.
    load(fullfile(blockFilePath, normBlockList{end}));
    eval(['block = ' normBlockList{end} ';']);
    clear(normBlockList{end});
    lastBlockSize = size(block{1}, 2);
    lastFrame = (blocks - 1) * blockSize + lastBlockSize - 1; 
    clear('block');
    
% Check the variables (the blocks are in a single file).
else
    
    % Load the worm information.
    load(wormFile, 'fps', 'lastFrame');
    
    % Check the variables.
    varName = 'fps';
    if ~exist(varName, 'var')
        error('worm2func:BadVariable',['Cannot load ''' varName ''' from ''' wormFile '''']);
    end
    varName = 'lastFrame';
    if ~exist(varName, 'var')
        error('worm2func:BadVariable', ['Cannot load ''' varName ''' from ''' wormFile '''']);
    end
end

% Correct the data types.
fps = double(fps);

% Compute the number of frames.
frames = lastFrame + 1;

%==========================================================================
%==========================================================================

%==========================================================================
%==========================================================================


%New inputs 
%------------------------------------------
fps = 20;



% Initialize the function state.
%
% Note: empirically I've found the values below achieve good signal.
%
% Furthermore ...
%
% Huang et al. in 2006, measure foraging frequencies for several worms and
% find the signal centered at roughly 4Hz. For N2 worms, they see a second
% signal at 10Hz but I find this value too close to the background noise
% present in segmentation. Visually inspecting the foraging signal, as the
% bend between the nose and neck, corroborates a roughly 4Hz signal. But,
% foraging usually encompasses only half to a quarter cycle. In other
% words, the worm bends it nose sharply and sometimes bends it back but a
% full wave, akin to a body bend, occurs far less frequently. Therefore I
% chose to measure angular speed for foraging.
%
% The body bend frequency is much easier to see. The N2 signal is clearly
% centered around 1/3Hz in both the literature and through visual
% inspection.
%
% I chose a high-frequency threshold of 4 frames. With 4 frames a 3-frame
% tick, resulting from segmentation noise, will be diluted by the
% additional frame.
%
% I chose a low-frequency threshold that requires at least half of the
% signal cycle to be present in the measurement window. In practice, this
% threshold appears to be unecessary as the data rarely, if ever, violates
% it.
minNoseWin     = round(0.1 * fps);
maxNoseWinTime = 15;
maxNoseWin     = round(maxNoseWinTime * fps);
noseState = struct( ...
    'minWin', minNoseWin, ...
    'maxWin', maxNoseWin, ...
    ... %'res', 2^11, ... % the FFT is quantized below this resolution
    'noseI', fliplr(1:4), ... % half the head (1/12 the worm)
    'neckI', fliplr(5:8)); % half the head (1/12 the worm)
    %'minFreq', 1 / (4 * maxNoseWinTime), ... % require at least 50% of the wave
    %'maxFreq', fps / 4, ... % with 4 frames we can resolve 75% of a wave
    %'minAmp', 15);
minBodyWinTime = .5;
minBodyWin = round(minBodyWinTime * fps);
maxBodyWinTime = 15;
maxBodyWin = round(maxBodyWinTime * fps);
bodyState = struct( ...
    'minWin', minBodyWin, ...
    'maxWin', maxBodyWin, ...
    'res',    2^14, ... % the FFT is quantized below this resolution
    'headI',  6:10, ... % centered at the head (1/6 the worm)
    'midI',   23:27, ... % centered at the middle of the worm
    'tailI',  40:44, ... % centered at the tail (1/6 the worm)
    'minFreq', 1 / (4 * maxBodyWinTime), ... % require at least 50% of the wave
    'maxFreq', fps / 4); % with 4 frames we can resolve 75% of a wave
    %'minAmp', 15); 
bendState = struct( ...
    'nose', noseState, ...
    'body', bodyState);

% Compute the bends.
win   = max(bendState.nose.maxWin / fps, bendState.body.maxWin / fps);

%???? - what is worm2func???
%I think it just handles partial blocks, ignore this for now ...
%Call function directly ...
%bendState - passed to calling function 
%
%
%   Also need data
%
% %
%     [dataInfo,blockInfo] = loadData(wormFile, fps, ...
%         startDataBlockFrame, endDataBlockFrame, backScaleFrames, ...
%         frontScaleFrames, isExtended, firstFrame, blocks, blockSize, ...
%         blockInfo);
%
% dataInfo = struct( ...
%     'startFrame',         startUseFrame, ...
%     'startDataFrame',     startFrame + backScale, ...
%     'startDataFrameI',    startFrameI, ...
%     'endFrame',           endUseFrame, ...
%     'endDataFrame',       endFrame - frontScale, ...
%     'endDataFrameI',      endFrameI, ...
%     'data', {data}, ...
%     'fps', fps);
%
%
%   ??? How are edge cases (1st and last frame) handled???
%
%   What is data?????
%

data  = worm2func(@h__bendFunc, bendState, wormFile, [], [], win, win);

%data would normally be a cell array with each 


bendData = [];
for i = 1:length(data)
    bendData = cat(1, bendData, data{i});
end

% Organize the data.
foragingBend = struct(...
    'amplitude', [], ...
    'angleSpeed', []);
bodyBend = struct(...
    'amplitude', [], ...
    'frequency', []);
bends = struct( ...
    'foraging', foragingBend, ...
    'head', bodyBend, ...
    'midbody', bodyBend, ...
    'tail', bodyBend);

% Clean up and offset the nose data.
%noseOff = round(bendState.nose.win * fps / 2);
%noseOff = round(bendState.nose.win * fps);
noseOff = 0;
noseAmps = [bendData.noseAmps];
noseFreqs = [bendData.noseFreqs];
if ventralMode > 1
    noseAmps = -noseAmps;
    noseFreqs = -noseFreqs;
end
bends.foraging.amplitude  = nan(1,frames);
bends.foraging.angleSpeed = nan(1,frames);
bends.foraging.amplitude((1 + noseOff):(length(noseAmps) + noseOff))   = noseAmps;
bends.foraging.angleSpeed((1 + noseOff):(length(noseFreqs) + noseOff)) = noseFreqs;

% Clean up and offset the head data.
%--------------------------------------------------------------------------
%bodyOff = round(bendState.body.win * fps / 2);
%bodyOff = round(bendState.body.win * fps);
bodyOff   = 0;
headAmps  = [bendData.headAmps];
headFreqs = [bendData.headFreqs];
bends.head.amplitude = nan(1,frames);
bends.head.frequency = nan(1,frames);
bends.head.amplitude((1 + bodyOff):(length(headAmps) + bodyOff))  = headAmps;
bends.head.frequency((1 + bodyOff):(length(headFreqs) + bodyOff)) = headFreqs;
if ~isempty(pausedMotion)
    bends.head.amplitude(pausedMotion) = NaN;
    bends.head.frequency(pausedMotion) = NaN;
end

% Clean up and offset the midbody data.
%--------------------------------------------------------------------------
midAmps  = [bendData.midAmps];
midFreqs = [bendData.midFreqs];
bends.midbody.amplitude = nan(1,frames);
bends.midbody.frequency = nan(1,frames);
bends.midbody.amplitude((1 + bodyOff):(length(midAmps) + bodyOff)) = ...
    midAmps;
bends.midbody.frequency((1 + bodyOff):(length(midFreqs) + bodyOff)) = ...
    midFreqs;
if ~isempty(pausedMotion)
    bends.midbody.amplitude(pausedMotion) = NaN;
    bends.midbody.frequency(pausedMotion) = NaN;
end

% Clean up and offset the tail data.
%--------------------------------------------------------------------------
tailAmps  = [bendData.tailAmps];
tailFreqs = [bendData.tailFreqs];
bends.tail.amplitude = nan(1,frames);
bends.tail.frequency = nan(1,frames);
bends.tail.amplitude((1 + bodyOff):(length(tailAmps) + bodyOff))  = tailAmps;
bends.tail.frequency((1 + bodyOff):(length(tailFreqs) + bodyOff)) = tailFreqs;
if ~isempty(pausedMotion)
    bends.tail.amplitude(pausedMotion) = NaN;
    bends.tail.frequency(pausedMotion) = NaN;
end
end



%% Compute the bend angles at the nose, head, midbody, and tail.
function bends = h__bendFunc(dataInfo, state)

%New inputs
%----------------------------
%fps
%data - 


%4 - skeletons
%5 - angles


%INPUTS
%==========================================================================
% Extract the nose and neck data.
noseSkeletons = data{4}(state.nose.noseI,:,:);
neckSkeletons = data{4}(state.nose.neckI,:,:);

% Extract the body data.
headBends = mean(data{5}(state.body.headI,:), 1);
midBends  = mean(data{5}(state.body.midI,:), 1);
tailBends = mean(data{5}(state.body.tailI,:), 1);
%==========================================================================

%start and end indices - may no longer be needed ...



% No worm data.
fps    = dataInfo.fps;
data   = dataInfo.data;
isData = data{1} == 's';
startI = dataInfo.startDataFrameI;
endI   = dataInfo.endDataFrameI;
if all(~isData)
    nanData = nan(1, endI - startI + 1);
    bends = struct( ...
        'noseAmps',     nanData, ...
        'noseFreqs',    nanData, ...
        'headAmps',     nanData, ...
        'headFreqs',    nanData, ...
        'midAmps',      nanData, ...
        'midFreqs',     nanData, ...
        'tailAmps',     nanData, ...
        'tailFreqs',    nanData);
    return;
end

% Find the start and end indices for missing data chunks.
isNotData        = ~isData;
isInterpNoseData = isNotData;
diffIsNotData    = diff(isNotData);
startNotDataI    = find(diffIsNotData == 1);
endNotDataI      = find(diffIsNotData == -1);

% Don't interpolate missing data at the very start and end.
if ~isempty(startNotDataI) && (isempty(endNotDataI) || startNotDataI(end) > endNotDataI(end))
    isInterpNoseData(startNotDataI(end):end) = false;
    startNotDataI(end) = [];
end
if ~isempty(endNotDataI) && ...
        (isempty(startNotDataI) || startNotDataI(1) > endNotDataI(1))
    isInterpNoseData(1:endNotDataI(1)) = false;
    endNotDataI(1) = [];
end

% Don't interpolate large missing chunks of data.
maxNoseInterp = 2 * state.nose.minWin - 1;
for i = 1:length(startNotDataI)
    if endNotDataI(i) - startNotDataI(i) > maxNoseInterp
        isInterpNoseData(startNotDataI(i):endNotDataI(i)) = false;
    end
end

%??? So it seems like data is 





% Interpolate the missing data.
%interpType = 'cubic';
%interpType = 'spline';
interpType  = 'linear';
dataI       = find(isData);
interpI     = find(isNotData);
noseInterpI = find(isInterpNoseData);
if ~isempty(interpI) && length(dataI) > 1
    
    % Interpolate the nose data.
    for i = 1:length(state.nose.noseI)
        noseSkeletons(i,1,noseInterpI) = interp1(dataI, squeeze(noseSkeletons(i,1,dataI)), noseInterpI, interpType, NaN);
        noseSkeletons(i,2,noseInterpI) = interp1(dataI, squeeze(noseSkeletons(i,2,dataI)), noseInterpI, interpType, NaN);
    end
    
    % Interpolate the neck data.
    for i = 1:length(state.nose.neckI)
        neckSkeletons(i,1,noseInterpI) = ...
            interp1(dataI, squeeze(neckSkeletons(i,1,dataI)), ...
            noseInterpI, interpType, NaN);
        neckSkeletons(i,2,noseInterpI) = ...
            interp1(dataI, squeeze(neckSkeletons(i,2,dataI)), ...
            noseInterpI, interpType, NaN);
    end
    
    % Interpolate the body data.
    headBends(interpI) = interp1(dataI, headBends(dataI), interpI, interpType, NaN);
    midBends(interpI)  = interp1(dataI, midBends(dataI), interpI, interpType, NaN);
    tailBends(interpI) = interp1(dataI, tailBends(dataI), interpI, interpType, NaN);
end

% Compute the nose bend angles.
noseDiffs = diff(noseSkeletons, 1, 1);
if size(noseDiffs, 1) > 1
    noseDiffs = mean(noseDiffs, 1);
end
noseAngles = squeeze(atan2(noseDiffs(:,2,:), noseDiffs(:,1,:)));
neckDiffs = diff(neckSkeletons, 1, 1);
if size(neckDiffs, 1) > 1
    neckDiffs = mean(neckDiffs, 1);
end
neckAngles = squeeze(atan2(neckDiffs(:,2,:), neckDiffs(:,1,:)));
noseBends  = (noseAngles - neckAngles)';
wrap       = noseBends > pi;
noseBends(wrap) = noseBends(wrap) - 2 * pi;
wrap       = noseBends < -pi;
noseBends(wrap) = noseBends(wrap) + 2 * pi;
noseBends  = noseBends * 180 / pi;




% Compute the worm bends.
[noseAmps,noseFreqs] = foragingData(noseBends, startI, endI, state.nose.minWin, fps);
[headAmps,headFreqs] = h__bendData(headBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps, state.body.res, false);
[midAmps,midFreqs] = h__bendData(midBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps, state.body.res, false);
[tailAmps,tailFreqs] = h__bendData(tailBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps,  state.body.res, false);

% Remove the noise.
%---------------------------------------------------------------
noseAmps(isnan(noseFreqs)) = NaN;
headFreqs(abs(headFreqs) > state.body.maxFreq) = NaN;
headFreqs(abs(headFreqs) < state.body.minFreq) = NaN;
headAmps(isnan(headFreqs)) = NaN;
midFreqs(abs(midFreqs) > state.body.maxFreq) = NaN;
midFreqs(abs(midFreqs) < state.body.minFreq) = NaN;
midAmps(isnan(midFreqs)) = NaN;
tailFreqs(abs(tailFreqs) > state.body.maxFreq) = NaN;
tailFreqs(abs(tailFreqs) < state.body.minFreq) = NaN;
tailAmps(isnan(tailFreqs)) = NaN;

% Organize the data.
bends = struct( ...
    'noseAmps',     noseAmps, ...
    'noseFreqs',    noseFreqs, ...
    'headAmps',     headAmps, ...
    'headFreqs',    headFreqs, ...
    'midAmps',      midAmps, ...
    'midFreqs',     midFreqs, ...
    'tailAmps',     tailAmps, ...
    'tailFreqs',    tailFreqs);


%bends
%.foraging
%   .amplitude
%   .frequency
%.head

end



%% Compute the bend amplitude and frequency.
function [amps,freqs] = h__bendData(data, startI, endI, minWinSize, maxWinSize, fps, fftRes, isCentering)
%
%   Inputs
%   =======================================================================
%   data   :
%   startI :
%   endI   :
%   minWinSize :
%   maxWinSize :
%   fps :
%   fftRes :
%   isCentering :
%   
%

% Compute the short-time Fourier transforms (STFT).
%dcThr = 0;
fftI  = 1:((fftRes + 1) / 2);
amps  = nan(1, endI - startI + 1);
freqs = nan(size(amps));
for i = 1:length(amps)
    
    % Pull out the time window.
    dataI = startI + i - 1;
    startDataWinI = max(1, dataI - maxWinSize);
    endDataWinI = min(length(data), dataI + maxWinSize);
    dataWin = data(startDataWinI:endDataWinI);
    dataWinI = dataI - startDataWinI + 1;
    
    % Center the signal.
    if isCentering
        dataWin = dataWin - mean(dataWin);
    end

    % Find the first zero crossing backwards.
    backZeroI = findZeroCrossing(dataWin, dataWinI, -1);
    if isempty(backZeroI)
        continue;
    end
    
    % Find the first zero crossing forwards.
    frontZeroI = findZeroCrossing(dataWin, dataWinI, 1);
    if isempty(frontZeroI)
        continue;
    end
    
    % Expand the zero-crossing window.
    %numZeros = 2;
    while frontZeroI - backZeroI + 1 < minWinSize %&& numZeros < 3
        if dataWinI - backZeroI < frontZeroI - dataWinI
            backZeroI = findZeroCrossing(dataWin, backZeroI - 1, -1);
            if isempty(backZeroI)
                break;
            end
        else
            frontZeroI = findZeroCrossing(dataWin, frontZeroI + 1, 1);
            if isempty(frontZeroI)
                break;
            end
        end
        %numZeros = numZeros + 1;
    end
    if isempty(backZeroI) || isempty(frontZeroI)
        continue;
    end
    
    % Center the window.
    %dataWinSize = round((frontZeroI - backZeroI - 1) / 2);
    dataWinSize = max(dataWinI - backZeroI, frontZeroI - dataWinI);
    backZeroI = dataWinI - dataWinSize;
    if backZeroI < 1
        continue;
    end
    frontZeroI = dataWinI + dataWinSize;
    if frontZeroI > length(dataWin)
        continue;
    end
    
    % Cut the window off at the zero crossings.
    dataWin = dataWin(backZeroI:frontZeroI);
    peakWinSize = round(sqrt(length(dataWin)));
    
    % Compute the real part of the STFT.
    fftData = fft(dataWin, fftRes);
    fftData = abs(fftData(fftI));
    
    % Find the peak frequency.
    [maxPeaks, maxPeaksI] = maxPeaksDist(fftData, peakWinSize);
    [~, maxI] = max(maxPeaks);
    maxPeakI = maxPeaksI(maxI);
    maxPeak = fftData(maxPeakI);
    
    % Find the peak bandwidth.
    [~, minPeaksI] = minPeaksDist(fftData, peakWinSize);
    peakStartI = minPeaksI(minPeaksI < maxPeakI);
    if ~isempty(peakStartI)
        peakStartI = peakStartI(end);
    end
    peakEndI = minPeaksI(minPeaksI > maxPeakI);
    if ~isempty(peakEndI)
        peakEndI = peakEndI(1);
    end
    
    % If the peak is not distinguished, we have no signal.
    peakBandThr = .5;
    peakEnergyThr = .5;
    if isempty(peakStartI) || isempty(peakEndI) || ...
            fftData(peakStartI) / maxPeak > peakBandThr || ...
            fftData(peakEndI) / maxPeak > peakBandThr || ...
            sum(fftData(peakStartI:peakEndI) .^ 2) / sum(fftData .^ 2) ...
            < peakEnergyThr
        maxPeakI = NaN;
    end
    
    % Convert the peak to a time frequency.
    if ~isnan(maxPeakI)
        dataSign = sign(mean(dataWin)); % sign the data
        amps(i) = (2 * fftData(maxPeakI) / length(dataWin)) * dataSign;
        freqs(i) = (fps / 2) * ((maxPeakI - 1) / (length(fftI) - 1)) * dataSign;
    end
    
    % Plot the STFT.
    % y = 2 * fftData / length(dataWin);
    % x = (fps / 2) * linspace(0, 1, length(fftI));
    % figure, plot(x, y);
    % xlabel('Frequency (Hz)');
    % ylabel('|Y(f)|');
    %
    % Plot the data signal.
    % figure, plot(dataWin);
end
end



%% Compute the foraging amplitude and angular speed.
function [amps speeds] = foragingData(data, startI, endI, minWinSize, fps)

% Initialize the amplitude and speed.
amps = nan(1, endI - startI + 1);
speeds = nan(size(amps));

% Clean up the signal with a gaussian filter.
if minWinSize > 0
    gaussFilter = gausswin(2 * minWinSize + 1) / minWinSize;
    data = conv(data, gaussFilter, 'same');
    data(1:minWinSize) = NaN;
    data((end - minWinSize + 1):end) = NaN;
end

% Compute the amplitudes between zero crossings.
dataSign = sign(data);
dataAmps = nan(1,length(data));
numAmps = 0;
for i = 1:(length(data) - 1)
    
    % Compute the amplitude for the region.
    % Note: data at the zero crossing has NaN (unknown) amplitude.
    if dataSign(i) ~= dataSign(i + 1);
        if dataSign(i) > 0
            dataAmps((i - numAmps):i) = max(data((i - numAmps):i));
        elseif dataSign(i) < 0
            dataAmps((i - numAmps):i) = min(data((i - numAmps):i));
        end
        
        % Reset the count.
        numAmps = 0;
        
    % Advance.
    else
        numAmps = numAmps + 1;
    end
end

% Compute the amplitude for the end region.
% Note: data at the zero crossing has NaN (unknown) amplitude.
if dataSign(end) > 0
    dataAmps((end - numAmps):end) = max(data((end - numAmps):end));
elseif dataSign(end) < 0
    dataAmps((end - numAmps):end) = min(data((end - numAmps):end));
end

% Compute the amplitude for our data.
amps = dataAmps(startI:endI);

% Compute the speed centered between the back and front foraging movements.
if startI == 1 && endI == length(data)
    dData = diff(data(startI:endI)) * fps;
    speeds(2:end-1) = (dData(1:(end - 1)) + dData(2:end)) / 2;
elseif startI == 1
    dData = diff(data(startI:(endI + 1))) * fps;
    speeds(2:end) = (dData(1:(end - 1)) + dData(2:end)) / 2;
elseif endI == length(data)
    dData = diff(data((startI - 1):endI)) * fps;
    speeds(1:(end-1)) = (dData(1:(end - 1)) + dData(2:end)) / 2;
else
    dData = diff(data((startI - 1):(endI + 1))) * fps;
    speeds(1:end) = (dData(1:(end - 1)) + dData(2:end)) / 2;
end
end



%% Compute the foraging amplitude and angular speed.
% function [amps speeds] = foragingData(data, startI, endI, ...
%     minWinSize, maxWinSize, fps, isCentering)
% 
% % Compute the amplitude and speed.
% amps = nan(1, endI - startI + 1);
% speeds = nan(size(amps));
% for i = 1:length(amps)
%     
%     % Pull out the time window.
%     dataI = startI + i - 1;
%     startDataWinI = max(1, dataI - maxWinSize);
%     endDataWinI = min(length(data), dataI + maxWinSize);
%     dataWin = data(startDataWinI:endDataWinI);
%     dataWinI = dataI - startDataWinI + 1;
%     
%     % The data point cannot be at the bounds of the time window.
%     if dataWinI <= 1 || dataWinI >= length(dataWin)
%         continue;
%     end
%     
%     % Center the signal.
%     if isCentering
%         dataWin = dataWin - mean(dataWin);
%     end
%     
%     % The data is at the zero crossing.
%     if dataWin(dataWinI) == 0
%         zeroI = dataWinI;
%         zeroMode = 0;
%         
%     % Find the nearest zero crossing.
%     else
%         backZeroI = findZeroCrossing(dataWin, dataWinI, -1);
%         frontZeroI = findZeroCrossing(dataWin, dataWinI, 1);
%         
%         % The zero crossing cannot be at the bounds of the time window.
%         if backZeroI == 1
%             backZeroI = [];
%         end
%         if frontZeroI == length(dataWin)
%             frontZeroI = [];
%         end
%         
%         % Is there a zero crossing in range?
%         if isempty(backZeroI) && isempty(frontZeroI)
%             continue;
%         end
%         
%         % Which zero crossing is nearest?
%         if isempty(backZeroI)
%             zeroI = frontZeroI;
%             zeroMode = 1;
%         elseif isempty(frontZeroI)
%             zeroI = backZeroI;
%             zeroMode = -1;
%         elseif dataWinI - backZeroI < frontZeroI - dataWinI
%             zeroI = backZeroI;
%             zeroMode = -1;
%         else
%             zeroI = frontZeroI;
%             zeroMode = 1;
%         end
%     end
%     
%     % Try the nearest zero crossing.
%     numTries = 2;
%     while numTries > 0
% 
%         % Find the peak on the back side of the zero crossing.
%         j = zeroI;
%         if zeroMode == -1
%             j = zeroI + 1;
%         end
%         minPeakI = findMinPeak(dataWin, j, minWinSize, -1);
%         maxPeakI = findMaxPeak(dataWin, j, minWinSize, -1);
%         if minPeakI >= zeroI
%             backPeakI = maxPeakI;
%         elseif maxPeakI >= zeroI
%             backPeakI = minPeakI;
%         else
%             backPeakI = min(minPeakI, maxPeakI);
%         end
%         
%         % Find the peak on the front side of the zero crossing.
%         j = zeroI;
%         if zeroMode == 1
%             j = zeroI - 1;
%         end
%         minPeakI = findMinPeak(dataWin, j, minWinSize, 1);
%         maxPeakI = findMaxPeak(dataWin, j, minWinSize, 1);
%         if minPeakI <= zeroI
%             frontPeakI = maxPeakI;
%         elseif maxPeakI <= zeroI
%             frontPeakI = minPeakI;
%         else
%             frontPeakI = min(minPeakI, maxPeakI);
%         end
%         
%         % Are we in the foraging window?
%         if dataWinI >= backPeakI && dataWinI <= frontPeakI
%             break;
%             
%         % Try again.
%         else
%             if numTries > 1
%                 
%                 % Try the front side of the zero crossing.
%                 numTries = numTries - 1;
%                 if zeroMode == -1 && ~isempty(frontZeroI)
%                     zeroI = frontZeroI;
%                     continue;
%                     
%                 % Try the back side of the zero crossing.
%                 elseif zeroMode == 1 && ~isempty(backZeroI)
%                     zeroI = backZeroI;
%                     continue;
%                 end
%             end
%             
%             % Fail.
%             numTries = 0;
%             zeroI = [];
%         end
%     end
%     
%     % Are we in a foraging window?
%     if isempty(zeroI)
%         continue;
%     end
%     
%     % Compute the foraging speed.
%     amps(i) = dataWin(frontPeakI) - dataWin(backPeakI);
%     speeds(i) = mean(diff(dataWin((dataWinI - 1):(dataWinI + 1)))) * fps;
% end
% end



%% Find the next zero crossing.
function zeroI = findZeroCrossing(data, i, increment)

% Is the index valid?
zeroI = [];
if i < 1 || i > length(data)
    return;
end

% Is the index at a zero crossing?
if data(i) == 0
    zeroI = i;
    return;
end

% Find the next zero crossing forward.
prevSign = sign(data(i));
i = i + increment;
if increment > 0
    while i <= length(data)

        % Did we cross zero?
        if sign(data(i)) ~= prevSign
            zeroI = i;
            return;
        end
        prevSign = sign(data(i));
        
        % Advance.
        i = i + increment;
    end
    
% Find the next zero crossing backward.
elseif increment < 0
    while i >= 1
        
        % Did we cross zero?
        if sign(data(i)) ~= prevSign
            zeroI = i;
            return;
        end
        prevSign = sign(data(i));
        
        % Advance.
        i = i + increment;
    end
end
end



%% Find the next minimum peak.
function peakI = findMinPeak(data, i, winSize, increment)

% Find the next minimum peak backward.
peakI = i;
i = i + increment;
if increment < 0
    while i >= 1 && peakI - i <= winSize
        
        % Record the new peak.
        if data(i) < data(peakI)
            peakI = i;
        end
        
        % Advance.
        i = i + increment;
    end
    
% Find the next minimum peak forward.
else
    while i <= length(data) && i - peakI <= winSize
        
        % Record the new peak.
        if data(i) < data(peakI)
            peakI = i;
        end
        
        % Advance.
        i = i + increment;
    end
end
end



%% Find the next maximum peak.
function peakI = findMaxPeak(data, i, winSize, increment)

% Find the next maximum peak backward.
peakI = i;
i = i + increment;
if increment < 0
    while i >= 1 && peakI - i <= winSize
        
        % Record the new peak.
        if data(i) > data(peakI)
            peakI = i;
        end
        
        % Advance.
        i = i + increment;
    end
    
% Find the next maximum peak forward.
else
    while i <= length(data) && i - peakI <= winSize
        
        % Record the new peak.
        if data(i) > data(peakI)
            peakI = i;
        end
        
        % Advance.
        i = i + increment;
    end
end
end
