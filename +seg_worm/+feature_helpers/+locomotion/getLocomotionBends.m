function bends = getLocomotionBends(sx,sy,angles,motion_mode,is_segmented_mask,ventral_mode)
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
%   Crawling. 
%   -----------------------
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
%
%   If the window between zero crossings is too small, the nearest zero
%   crossing is assumed to be noise and we search for the next available
%   zero crossing in its respective direction. If the window is too big,
%   crawling is marked undefined at the frame. 
%
%   Once an appropriate window has been found, the window is extended in
%   order to center the frame and measure instantaneous crawling by
%   ensuring that the distance on either side to respective zero crossings
%   is identical. If the distances are not identical, the distance of the
%   larger side is used in place of the zero-crossing distance of the
%   smaller side in order to expand the small side and achieve a symmetric
%   window, centered at the frame of interest.
% 
%   We use a Fourier transform to measure the amplitude and frequency
%   within the window described above. The largest peak within the
%   transform is chosen for the crawling amplitude and frequency. If the
%   troughs on either side of the peak exceed 1/2 its height, the peak is
%   rejected for being unclear and crawling is marked as undefined at the
%   frame. Similarly, if the integral between the troughs is less than half
%   the total integral, the peak is rejected for being weak.
%




%New inputs 
%------------------------------------------
FPS = 20;

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


minNoseWin     = round(0.1 * FPS);
maxNoseWinTime = 15;
maxNoseWin     = round(maxNoseWinTime * FPS);
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
minBodyWin = round(minBodyWinTime * FPS);
maxBodyWinTime = 15;
maxBodyWin = round(maxBodyWinTime * FPS);
bodyState = struct( ...
    'minWin', minBodyWin, ...
    'maxWin', maxBodyWin, ...
    'res',    2^14, ... % the FFT is quantized below this resolution
    'headI',  6:10, ... % centered at the head (1/6 the worm)
    'midI',   23:27, ... % centered at the middle of the worm
    'tailI',  40:44, ... % centered at the tail (1/6 the worm)
    'minFreq', 1 / (4 * maxBodyWinTime), ... % require at least 50% of the wave
    'maxFreq', FPS / 4); % with 4 frames we can resolve 75% of a wave
    %'minAmp', 15); 
bendState = struct( ...
    'nose', noseState, ...
    'body', bodyState);




%{
    
% Compute the bends.
win   = max(bendState.nose.maxWin / FPS, bendState.body.maxWin / FPS);    
    
wormFile = 'F:\worm_data\segworm_data\video\testing_with_GUI\.data\mec-4 (u253) off food x_2010_04_21__17_19_20__1_seg\normalized\segNormInfo.mat';

%This took for ever (30 - 40 s?), don't need to rerun, saved results,
%unless I am debugging the object form ...
data = seg_worm.w.util.worm2func(@h__bendFuncOld, bendState, wormFile, [], [], win, win);

% old_data = data;
% save('wtf','old_data')
%}

data = h__bendFunc(bendState,FPS,angles,is_segmented_mask);

keyboard



%data would normally be a cell array with each 
%element containing a chunk or block of data from a subset of frames ...

%Turn the cell array into a structure array
bendData = [];
for i = 1:length(data)
    bendData = cat(1, bendData, data{i});
end





%The following all looks very similar ...
%--------------------------------------------------------------------------

%TODO: Make below into a function ...

%Actually, depending upon the format, all this really does 
%is 

% Clean up and offset the head data.
%--------------------------------------------------------------------------
%bodyOff = round(bendState.body.win * fps / 2);
%bodyOff = round(bendState.body.win * fps);

bodyOff   = 0; %This seems like we 



end

function h__getCrawlingStructs()


bodyBend = struct(...
    'amplitude', [], ...
    'frequency', []);
bends = struct( ...
    'foraging', foragingBend, ...
    'head',     bodyBend, ...
    'midbody',  bodyBend, ...
    'tail',     bodyBend);

%JAH: Unfinished ...

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
bends.midbody.amplitude((1 + bodyOff):(length(midAmps) + bodyOff))  = midAmps;
bends.midbody.frequency((1 + bodyOff):(length(midFreqs) + bodyOff)) = midFreqs;

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


function bends = h__bendFunc(state,fps,angles,is_segmented_mask)
%
%
%   bends = h__bendFunc(dataInfo, state)
%New inputs
%----------------------------
%fps
%data - 


%OLD
% % % dataInfo:
% % %          startFrame: 0
% % %      startDataFrame: 0
% % %     startDataFrameI: 1
% % %            endFrame: 925
% % %        endDataFrame: 499
% % %       endDataFrameI: 500
% % %                data: {12x1 cell}
% % %                 fps: 25.8398
%




%==========================================================================
% state.nose
%     minWin: 2
%     maxWin: 300
%      noseI: [4 3 2 1]
%      neckI: [8 7 6 5]
% state.body
%      minWin: 10
%      maxWin: 300
%         res: 16384
%       headI: [6 7 8 9 10]
%        midI: [23 24 25 26 27]
%       tailI: [40 41 42 43 44]
%     minFreq: 0.0167
%     maxFreq: 5

%4 - skeletons
%5 - angles


%CONSTANTS
%==========================================================================

SI = seg_worm.skeleton_indices;

%headI   6:10
%midI    23:27
%tailI   40:44

%INPUTS
%==========================================================================

% Extract the body data.
headBends = mean(angles(state.body.headI,:), 1);
midBends  = mean(angles(state.body.midI,:), 1);
tailBends = mean(angles(state.body.tailI,:), 1);
%==========================================================================

%start and end indices - may no longer be needed ...

% No worm data.

isData = is_segmented_mask;
startI = 1;
endI   = length(is_segmented_mask);



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

%Interpolation ....
%==========================================================================
% Find the start and end indices for missing data chunks.
isNotData        = ~isData;



%TODO: Do one set of commands, three times, i.e. remove the set of 
%3 for the calls ...

%bends = headBends, midBends, tailBends
%
%   Then run code and get amps and freqs ...


% Interpolate the missing data.
interpType  = 'linear';
dataI       = find(isData);
interpI     = find(isNotData);
if ~isempty(interpI) && length(dataI) > 1
    headBends(interpI) = interp1(dataI, headBends(dataI), interpI, interpType, NaN);
    midBends(interpI)  = interp1(dataI, midBends(dataI), interpI, interpType, NaN);
    tailBends(interpI) = interp1(dataI, tailBends(dataI), interpI, interpType, NaN);
end
%==========================================================================


%JAH: At this point ...

keyboard

[headAmps,headFreqs] = h__bendData(headBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps, state.body.res, false);
[midAmps,midFreqs] = h__bendData(midBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps, state.body.res, false);
[tailAmps,tailFreqs] = h__bendData(tailBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps,  state.body.res, false);

% Remove the noise.
%---------------------------------------------------------------


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
    'noseFreqs',    noseFreqs, ...
    'headAmps',     headAmps, ...
    'headFreqs',    headFreqs, ...
    'midAmps',      midAmps, ...
    'midFreqs',     midFreqs, ...
    'tailAmps',     tailAmps, ...
    'tailFreqs',    tailFreqs);

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
    endDataWinI   = min(length(data), dataI + maxWinSize);
    dataWin       = data(startDataWinI:endDataWinI);
    dataWinI      = dataI - startDataWinI + 1;
    
    % Center the signal.
    if isCentering
        dataWin = dataWin - mean(dataWin);
    end

    % Find the first zero crossing backwards.
    backZeroI = h__findZeroCrossing(dataWin, dataWinI, -1);
    if isempty(backZeroI)
        continue;
    end
    
    % Find the first zero crossing forwards.
    frontZeroI = h__findZeroCrossing(dataWin, dataWinI, 1);
    if isempty(frontZeroI)
        continue;
    end
    
    % Expand the zero-crossing window.
    %numZeros = 2;
    while frontZeroI - backZeroI + 1 < minWinSize %&& numZeros < 3
        if dataWinI - backZeroI < frontZeroI - dataWinI
            backZeroI = h__findZeroCrossing(dataWin, backZeroI - 1, -1);
            if isempty(backZeroI)
                break;
            end
        else
            frontZeroI = h__findZeroCrossing(dataWin, frontZeroI + 1, 1);
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
function [amps,speeds] = h__foragingData(data, startI, endI, minWinSize, fps)

% Initialize the amplitude and speed.
amps = nan(1, endI - startI + 1);
speeds = nan(size(amps));

% Clean up the signal with a gaussian filter.
if minWinSize > 0
    gaussFilter = gausswin(2 * minWinSize + 1) / minWinSize;
    data        = conv(data, gaussFilter, 'same');
    data(1:minWinSize) = NaN;
    data((end - minWinSize + 1):end) = NaN;
end

% Compute the amplitudes between zero crossings.
dataSign = sign(data);
dataAmps = nan(1,length(data));
numAmps  = 0;
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










%% Find the next zero crossing.
function zeroI = h__findZeroCrossing(data, i, increment)
%
%
%   Inputs
%   =======================================================================
%   data :
%   i    :
%   increment :

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





%Start of the Old Code

%{





%==========================================================================
%==========================================================================
%==========================================================================

%% Compute the foraging amplitude and angular speed.
function [amps,speeds] = foragingData(data, startI, endI, minWinSize, fps)

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

%% Compute the bend amplitude and frequency.
function [amps,freqs] = h__bendDataOld(data, startI, endI, minWinSize, maxWinSize, fps, fftRes, isCentering)
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
    endDataWinI   = min(length(data), dataI + maxWinSize);
    dataWin       = data(startDataWinI:endDataWinI);
    dataWinI      = dataI - startDataWinI + 1;
    
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
    [maxPeaks, maxPeaksI] = seg_worm.util.maxPeaksDist(fftData, peakWinSize,true,-Inf);
    [~, maxI] = max(maxPeaks);
    maxPeakI  = maxPeaksI(maxI);
    maxPeak   = fftData(maxPeakI);
    
    
    
    
    % Find the peak bandwidth.
    [~, minPeaksI] = seg_worm.util.maxPeaksDist(fftData, peakWinSize,false,Inf);
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

%% Compute the bend angles at the nose, head, midbody, and tail.
function [bends,state] = h__bendFuncOld(dataInfo, state)

% No worm data.
fps = dataInfo.fps;
data = dataInfo.data;
isData = data{1} == 's';
startI = dataInfo.startDataFrameI;
endI = dataInfo.endDataFrameI;
if all(~isData)
    nanData = nan(1, endI - startI + 1);
    bends = struct( ...
        'noseAmps', nanData, ...
        'noseFreqs', nanData, ...
        'headAmps', nanData, ...
        'headFreqs', nanData, ...
        'midAmps', nanData, ...
        'midFreqs', nanData, ...
        'tailAmps', nanData, ...
        'tailFreqs', nanData);
    return;
end

% Find the start and end indices for missing data chunks.
isNotData = ~isData;
isInterpNoseData = isNotData;
diffIsNotData = diff(isNotData);
startNotDataI = find(diffIsNotData == 1);
endNotDataI = find(diffIsNotData == -1);

% Don't interpolate missing data at the very start and end.
if ~isempty(startNotDataI) && ...
        (isempty(endNotDataI) || startNotDataI(end) > endNotDataI(end))
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

% Extract the nose and neck data.
noseSkeletons = data{4}(state.nose.noseI,:,:);
neckSkeletons = data{4}(state.nose.neckI,:,:);

% Extract the body data.
headBends = mean(data{5}(state.body.headI,:), 1);
midBends  = mean(data{5}(state.body.midI,:), 1);
tailBends = mean(data{5}(state.body.tailI,:), 1);

% Interpolate the missing data.
%interpType = 'cubic';
%interpType = 'spline';
interpType = 'linear';
dataI = find(isData);
interpI = find(isNotData);
noseInterpI = find(isInterpNoseData);
if ~isempty(interpI) && length(dataI) > 1
    
    % Interpolate the nose data.
    for i = 1:length(state.nose.noseI)
        noseSkeletons(i,1,noseInterpI) = ...
            interp1(dataI, squeeze(noseSkeletons(i,1,dataI)), ...
            noseInterpI, interpType, NaN);
        noseSkeletons(i,2,noseInterpI) = ...
            interp1(dataI, squeeze(noseSkeletons(i,2,dataI)), ...
            noseInterpI, interpType, NaN);
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
    headBends(interpI) = ...
        interp1(dataI, headBends(dataI), interpI, interpType, NaN);
    midBends(interpI) = ...
        interp1(dataI, midBends(dataI), interpI, interpType, NaN);
    tailBends(interpI) = ...
        interp1(dataI, tailBends(dataI), interpI, interpType, NaN);
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
noseBends = (noseAngles - neckAngles)';
wrap = noseBends > pi;
noseBends(wrap) = noseBends(wrap) - 2 * pi;
wrap = noseBends < -pi;
noseBends(wrap) = noseBends(wrap) + 2 * pi;
noseBends = noseBends * 180 / pi;

% Compute the worm bends.
[noseAmps noseFreqs] = foragingData(noseBends, startI, endI, ...
    state.nose.minWin, fps);
[headAmps headFreqs] = h__bendDataOld(headBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps, state.body.res, false);
[midAmps midFreqs] = h__bendDataOld(midBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps, state.body.res, false);
[tailAmps tailFreqs] = h__bendDataOld(tailBends, startI, endI, ...
    state.body.minWin, state.body.maxWin, fps,  state.body.res, false);

% Remove the noise.
%noseFreqs(noseFreqs > state.nose.maxThr) = NaN;
%noseFreqs(noseFreqs < state.nose.minThr) = NaN;
noseAmps(isnan(noseFreqs)) = NaN;
headFreqs(abs(headFreqs) > state.body.maxFreq) = NaN;
headFreqs(abs(headFreqs) < state.body.minFreq) = NaN;
%headFreqs(abs(headAmps) < state.body.minAmp) = NaN;
headAmps(isnan(headFreqs)) = NaN;
midFreqs(abs(midFreqs) > state.body.maxFreq) = NaN;
midFreqs(abs(midFreqs) < state.body.minFreq) = NaN;
%midFreqs(abs(midAmps) < state.body.minAmp) = NaN;
midAmps(isnan(midFreqs)) = NaN;
tailFreqs(abs(tailFreqs) > state.body.maxFreq) = NaN;
tailFreqs(abs(tailFreqs) < state.body.minFreq) = NaN;
%tailFreqs(abs(tailAmps) < state.body.minAmp) = NaN;
tailAmps(isnan(tailFreqs)) = NaN;

% Organize the data.
bends = struct( ...
    'noseAmps', noseAmps, ...
    'noseFreqs', noseFreqs, ...
    'headAmps', headAmps, ...
    'headFreqs', headFreqs, ...
    'midAmps', midAmps, ...
    'midFreqs', midFreqs, ...
    'tailAmps', tailAmps, ...
    'tailFreqs', tailFreqs);
end



%}