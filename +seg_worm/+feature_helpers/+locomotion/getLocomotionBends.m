function bends = getLocomotionBends(angles,is_paused,is_segmented_mask)
%WORMBENDS Compute the temporal bending frequency at the nose, head,
%midbody, and tail.
%
%
%   bends =
%   seg_worm.feature_helpers.locomotion.getLocomotionBends(angles,is_paused,is_segmented_mask)
%
%
%   Old Name: wormBends.m
%
%   Inputs
%   =======================================================================
%   angles            :
%   is_paused         :
%   is_segmented_mask :
%
%
%   Outputs:
%   =======================================================================
%   bends
%       .head
%           .amplitude
%           .frequency
%       .mid
%       .tail
%
%   Nature Methods Description - Crawling and Foraging (might split)
%   =======================================================================
%   Crawling.
%   -----------------------
%   Worm crawling is expressed as both an amplitude and frequency
%   (Supplementary Fig. 4e). We measure these features instantaneously at
%   the head, midbody, and tail. The amplitude and frequency are signed
%   negatively whenever the worm�s ventral side is contained within the
%   concave portion of its instantaneous bend.
%
%   Crawling is only measured during forward and backward motion states.
%   The worm bend mean angles (described in the section on �Posture�) show
%   a roughly periodic signal as the crawling wave travels along the worm�s
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
%   bend angle mean � the location where the measured body part (head,
%   midbody, or tail) must have hit a flat posture (a supplementary bend
%   angle of 0�). This guarantees that we are observing half a cycle for
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
%   TODO:
%   - finish documentation
%   - handle FPS input



%New inputs
%------------------------------------------
FPS = 25;

% Initialize the function state.
%
% Note: empirically I've found the values below achieve good signal.
%
% Furthermore ...
%
% The body bend frequency is much easier to see (than foraging). The N2
% signal is clearly centered around 1/3Hz in both the literature and
% through visual inspection.
%
% I chose a high-frequency threshold of 4 frames. With 4 frames a 3-frame
% tick, resulting from segmentation noise, will be diluted by the
% additional frame.
%
% I chose a low-frequency threshold that requires at least half of the
% signal cycle to be present in the measurement window. In practice, this
% threshold appears to be unecessary as the data rarely, if ever, violates
% it.


%Setup Options
%------------------------------------------------------------------
minBodyWinTime = .5;
minBodyWin     = round(minBodyWinTime * FPS);
maxBodyWinTime = 15;
maxBodyWin     = round(maxBodyWinTime * FPS);
options = struct( ...
    'minWin',   minBodyWin, ...
    'maxWin',   maxBodyWin, ...
    'res',      2^14, ... % the FFT is quantized below this resolution
    'headI',    6:10, ... % centered at the head (1/6 the worm)
    'midI',     23:27, ... % centered at the middle of the worm
    'tailI',    40:44, ... % centered at the tail (1/6 the worm)
    'minFreq',  1 / (4 * maxBodyWinTime), ... % require at least 50% of the wave
    'maxFreq',  FPS / 4, ... % with 4 frames we can resolve 75% of a wave
    'peakBandThr',0.5,...
    'peakEnergyThr',0.5);
%Why are the indices used that are used ????
%NOTE: These indices are not symettric, 
%Should we use the skeleton indices instead, or move these to the skeleton
%indices????
%SI = seg_worm.skeleton_indices;

% No worm data.
%------------------------------------
if ~any(is_segmented_mask)
    nanData = nan(1, length(is_segmented_mask));
    bend_struct = struct('frequency',nanData,'amplitude',nanData);
    bends = struct( ...
        'head',     bend_struct, ...
        'mid',      bend_struct, ...
        'tail',     bend_struct);
    return
end

%Set things up for the loop
%------------------------------------
section     = {'head' 'mid' 'tail'};
all_indices = {options.headI options.midI options.tailI};

%Initialize interpolation - we'll interpolate over all frames but no extrapolation
%--------------------------------------------------------------------------
interpType  = 'linear';
dataI       = find(is_segmented_mask);
interpI     = find(~is_segmented_mask);

bends = struct();
for iBend = 1:3
    
    cur_indices     = all_indices{iBend};
    avg_bend_angles = mean(angles(cur_indices,:), 1);

    if ~isempty(interpI) && length(dataI) > 1
        avg_bend_angles(interpI) = interp1(dataI, avg_bend_angles(dataI), interpI, interpType, NaN);
    end
    
    [amps,freqs] = h__getBendData(avg_bend_angles,FPS,options,is_paused);
    
    bends.(section{iBend}).frequency = freqs;
    bends.(section{iBend}).amplitude = amps;
end
    
end



function [amps,freqs] = h__getBendData(avg_bend_angles, fps, options, is_paused)
%h__getBendData    Compute the bend amplitude and frequency.
%
%   h__getBendData(avg_bend_angles, fps, options, is_paused)


INIT_MAX_I_FOR_BANDWIDTH = 2000; %TODO: Relate this to a real frequency ...
%and pass it in from higher up. This number is not important to the
%outcome, only to the speed in which the function runs. We try and find the
%bandwidth within this number of samples. 


%Options extraction
%----------------------------------
min_window = options.minWin;
max_window = options.maxWin;
max_freq   = options.maxFreq;
min_freq   = options.minFreq;
fft_n_samples = options.res;
peakBandThr   = options.peakBandThr;
peakEnergyThr = options.peakEnergyThr;

%TODO: This needs to be cleaned up ...
[back_zeros_I,front_zeros_I] = h__getBoundingZeroIndices(avg_bend_angles,min_window);

n_frames = length(avg_bend_angles);

left_distances  = (1:n_frames) - back_zeros_I;
right_distances = front_zeros_I - (1:n_frames);
half_distances  = max(left_distances,right_distances);


left_bounds  = (1:n_frames) - half_distances;
right_bounds = (1:n_frames) + half_distances;

%Compute conditions by which we will ignore frames:
%--------------------------------------------------------------------------
%- frame is not bounded on both sides by a sign change
%- avg_bend_angles is NaN, this will only happen on the edges because we
%   interpolate over the other frames ... (we just don't extrapolate)
%- the sign change region is too large
%- the bounds we settle on exceed the data region
%- mode segmentation determined the frame was a paused frame
%??? - what about large NaN regions, are those paused regions???

is_bad_mask  = back_zeros_I == 0 | front_zeros_I == 0 | isnan(avg_bend_angles) | half_distances > max_window;
is_bad_mask  = is_bad_mask | left_bounds < 1 | right_bounds > n_frames | is_paused;


% Compute the short-time Fourier transforms (STFT).
fft_max_I   = fft_n_samples / 2; %Maximum index to keep for frequency analysis
freq_scalar = (fps / 2) * 1/(fft_max_I - 1);

n_frames = length(avg_bend_angles);

amps  = NaN(1,n_frames);
freqs = NaN(1,n_frames);
for iFrame = find(~is_bad_mask)
    
    dataWin = avg_bend_angles(left_bounds(iFrame):right_bounds(iFrame));
    data_win_length = length(dataWin);
    
    %fft frequency and bandwidth
    %----------------------------------------------------------------------
    % Compute the real part of the STFT.
    %These two steps take a lot of time ...
    fftData = fft(dataWin, fft_n_samples);
    fftData = abs(fftData(1:fft_max_I));
    
    % Find the peak frequency.
    [maxPeak,maxPeakI] = max(fftData);
            
    %NOTE: If this is true, we'll never bound the peak on the left ...
    if maxPeakI == 1
        continue
    end
    
    %TODO: Not sure if this value is correct ...
    unsigned_freq = freq_scalar*(maxPeakI - 1);
    
    if unsigned_freq > max_freq || unsigned_freq < min_freq
       continue 
    end
    
    [peakStartI,peakEndI] = h__getBandwidth(data_win_length,fftData,maxPeakI,INIT_MAX_I_FOR_BANDWIDTH);

    %Store data
    %----------------------------------------------------------------------
    if ~(   isempty(peakStartI)                           || ...
            isempty(peakEndI)                             || ...
            fftData(peakStartI)/maxPeak  > peakBandThr    || ...
            fftData(peakEndI)/maxPeak   > peakBandThr     || ...
            sum(fftData(peakStartI:peakEndI) .^ 2) / sum(fftData .^ 2) < peakEnergyThr)

        % Convert the peak to a time frequency.
        dataSign      = sign(mean(dataWin)); % sign the data
        amps(iFrame)  = (2 * fftData(maxPeakI) / data_win_length) * dataSign;
        freqs(iFrame) = unsigned_freq * dataSign;
    end
    
end

end

function [peakStartI,peakEndI] = h__getBandwidth(data_win_length,fftData,maxPeakI,INIT_MAX_I_FOR_BANDWIDTH)

peakWinSize = round(sqrt(data_win_length));

    % Find the peak bandwidth.
    %----------------------------------------------------------------------
    %The goal is to mind minimum 'peaks' that border the frequency
    %response.
    [~, minPeaksI] = seg_worm.util.maxPeaksDist(fftData(1:INIT_MAX_I_FOR_BANDWIDTH), peakWinSize,false,Inf);
    
    peakStartI = minPeaksI(find(minPeaksI < maxPeakI,1));
    peakEndI   = minPeaksI(find(minPeaksI > maxPeakI,1));
    
    if isempty(peakEndI) || peakEndI + peakWinSize >= INIT_MAX_I_FOR_BANDWIDTH   
        [~, minPeaksI] = seg_worm.util.maxPeaksDist(fftData, peakWinSize,false,Inf);

        peakStartI = minPeaksI(find(minPeaksI < maxPeakI,1));
        peakEndI   = minPeaksI(find(minPeaksI > maxPeakI,1));      
    end


end

function [amps,freqs] = h__temp(avg_bend_angles, minWinSize, maxWinSize, fps, fftRes,max_freq,min_freq)


INIT_MAX_I_FOR_BANDWIDTH = 2000; %TODO: Relate this to a real frequency ...
%and pass it in from higher up ...

peakBandThr   = .5;
peakEnergyThr = .5;

% Compute the short-time Fourier transforms (STFT).
%dcThr = 0;

fft_max_I   = fftRes / 2; %Maximum index to keep for frequency analysis
freq_scalar = (fps / 2) * 1/(fft_max_I - 1);

n_frames = length(avg_bend_angles);

amps  = NaN(1,n_frames);
freqs = NaN(1,n_frames);
for iFrame = 1:n_frames
    
    % Pull out the time window.
    %---------------------------------------------------
    startDataWinI = max(1, iFrame - maxWinSize);
    endDataWinI   = min(length(avg_bend_angles), iFrame + maxWinSize);
    
    
    %I'd like to delay this step until computing the fft ...
    dataWin       = avg_bend_angles(startDataWinI:endDataWinI);
    
    %?????
    dataWinI      = iFrame - startDataWinI + 1;
    
    %----------------------------------------------------------------------

    [backZeroI,frontZeroI] = h__getZeroCrossings(dataWin,dataWinI,minWinSize);
    if isempty(backZeroI) || isempty(frontZeroI)
        continue;
    end
        
    %----------------------------------------------------------------------
    

    % Center the window.
    %----------------------------------------------------------------------
    dataWinSize = max(dataWinI - backZeroI, frontZeroI - dataWinI);
    
    backZeroI   = dataWinI - dataWinSize;
    if backZeroI < 1
        continue;
    end
    
    frontZeroI = dataWinI + dataWinSize;
    if frontZeroI > length(dataWin)
        continue;
    end
    
    % Cut the window off at the zero crossings.
    dataWin     = dataWin(backZeroI:frontZeroI);
    peakWinSize = round(sqrt(length(dataWin)));
    
    
    %fft frequency and bandwidth
    %----------------------------------------------------------------------
    % Compute the real part of the STFT.
    %These two steps take a lot of time ...
    fftData = fft(dataWin, fftRes);
    fftData = abs(fftData(1:fft_max_I));
    
    % Find the peak frequency.
    [maxPeak,maxPeakI] = max(fftData);
            
    %NOTE: If this is true, we'll never bound the peak on the left ...
    if maxPeakI == 1
        continue
    end
    
    %TODO: Not sure if this value is correct ...
    unsigned_freq = freq_scalar*(maxPeakI - 1);
    
    if unsigned_freq > max_freq || unsigned_freq < min_freq
       continue 
    end
    
    % Find the peak bandwidth.
    %----------------------------------------------------------------------
    %The goal is to mind minimum 'peaks' that border the frequency
    %response.
    [~, minPeaksI] = seg_worm.util.maxPeaksDist(fftData(1:INIT_MAX_I_FOR_BANDWIDTH), peakWinSize,false,Inf);
    
    peakStartI = minPeaksI(find(minPeaksI < maxPeakI,1));
    peakEndI   = minPeaksI(find(minPeaksI > maxPeakI,1));
    
    if isempty(peakEndI) || peakEndI + peakWinSize >= INIT_MAX_I_FOR_BANDWIDTH   
        [~, minPeaksI] = seg_worm.util.maxPeaksDist(fftData, peakWinSize,false,Inf);

        peakStartI = minPeaksI(find(minPeaksI < maxPeakI,1));
        peakEndI   = minPeaksI(find(minPeaksI > maxPeakI,1));      
    end
    
    %Store data
    %----------------------------------------------------------------------
    if ~(   isempty(peakStartI)                           || ...
            isempty(peakEndI)                             || ...
            fftData(peakStartI)/maxPeak  > peakBandThr    || ...
            fftData(peakEndI)/maxPeak   > peakBandThr     || ...
            sum(fftData(peakStartI:peakEndI) .^ 2) / sum(fftData .^ 2) < peakEnergyThr)

        % Convert the peak to a time frequency.
        dataSign      = sign(mean(dataWin)); % sign the data
        amps(iFrame)  = (2 * fftData(maxPeakI) / length(dataWin)) * dataSign;
        freqs(iFrame) = unsigned_freq * dataSign;
    end
    
end

end


function [back_zeros_I,front_zeros_I] = h__getBoundingZeroIndices(avg_bend_angles,minWinSize)
%
%
%
%   The goal of this function is to bound each index by 

%TODO: Finish documentation

sign_change_I  = find(sign(avg_bend_angles(2:end)) ~= sign(avg_bend_angles(1:end-1)));
n_sign_changes = length(sign_change_I);

%backward - at sign changes - don't subtract or add
%forward  - we need to add 1

%???? How to propagate indices indices????

%Let's say we have indices 3  6  9
%What we need ...
%        1 2 3 4 5 6 7 9 10
%Left  = 0 0 0 3 3 3 6 6 6  %At 4, we have a zero at 3 to the left ...
%Right = 2 2 5 5 5 8 8 8

% dSignChange = diff(sign_change_I);

n_frames = length(avg_bend_angles);

%For each element in the array, these values indicate which sign change
%index to use ...
left_sign_change_I  = zeros(1,n_frames);
left_sign_change_I(sign_change_I + 1) = 1;
left_sign_change_I = cumsum(left_sign_change_I);

right_sign_change_I    = zeros(1,n_frames);
right_sign_change_I(sign_change_I(1:end-1)+1) = 1;
right_sign_change_I(1) = 1;
right_sign_change_I    = cumsum(right_sign_change_I);
right_sign_change_I(sign_change_I(end)+1:end) = 0; %Nothing to the right of the last change ...

%These are the actual indices that each sign crossing index points to
%
%We keep these separate as it makes it easier to go the next value, by
%incrementing the pointer index, rather than doing a search
left_values  = sign_change_I;
right_values = sign_change_I + 1; 


back_zeros_I  = zeros(1,n_frames);
front_zeros_I = zeros(1,n_frames);

for iFrame = 1:n_frames
    
    cur_left_index = left_sign_change_I(iFrame);
    cur_right_index = right_sign_change_I(iFrame);
    if left_sign_change_I(iFrame) == 0 || right_sign_change_I(iFrame) == 0
        continue
    end
    
    back_zero_I  = left_values(cur_left_index);
    front_zero_I = right_values(cur_right_index);
    
    use_values = true;
    % Expand the zero-crossing window.
    while front_zero_I - back_zero_I + 1 < minWinSize
        
        %If the distance of the zero crossing going forward is further
        %from the current index than the backward crosssing is, then expand
        %the backward crossing
        if iFrame - back_zero_I < front_zero_I - iFrame
            %Then expand backwards
            cur_left_index = cur_left_index - 1;
            if cur_left_index == 0
                use_values = false;
                break
            end
            back_zero_I  = left_values(cur_left_index);
        else
            %Expand forwards 
            cur_right_index = cur_right_index + 1;
            if cur_right_index > n_sign_changes
                use_values = false;
                break
            end
            front_zero_I = right_values(cur_right_index);
        end
    end
    
    if use_values
        back_zeros_I(iFrame)  = back_zero_I;
        front_zeros_I(iFrame) = front_zero_I;
    end
end

start_I = back_zeros_I;
stop_I  = front_zeros_I;


end