function numKinks = wormKinks(worm_angles)
%WORMKINKS Compute the kinks in a worm.
%
%
%   num_kinks = seg_worm.feature_helpers.posture.wormKinks(worm_angles)
%
%
%   [NUMKINKS KINKANGLES KINKINDICES] = WORMKINKS(WORMANGLES, LENGTHTHR,
%                                                 AMPTHR, ISSMOOTHING,
%                                                 WORMSKELETONS)
%
%   Inputs:
%       wormAngles    - the worm(s) bend angles at each skeleton point
%       lengthThr     - the bend segment length threshold, shorter segments
%                       are considered noise, not bends;
%                       the default is 1/12 of the worm
%       ampThr        - the bend amplitude threshold, smaller amplitudes
%                       are considered noise, not bends;
%                       the default is 0 degrees
%       issmoothing   - are we smoothing the worm angles? If so, the angles
%                       are convolved with a gaussian filter of the same
%                       size as the length threshold; the default is true
%       wormSkeletons - a 3-D matrix of worms (the first 2 dimensions are
%                       the x and y coordinates and the 3rd dimension is
%                       worms); when present, the worm kinks are displayed
%
%   Outputs:
%       numKinks - the number of kinks
%       indices  - the indices of the kinks
%
%   Nature Methods Description
%   =======================================================================
%   Bend Count - kinks
%   ---------------------------
%   The bend count is a rough measure of the number of bends along the
%   worm. The supplementary skeleton angles are measured during
%   segmentation and signed to reflect their dorsal-ventral orientation.
%
%   These angles are convolved with a Gaussian filter, 1/12 the length of
%   the skeleton, with a width defined by the Matlab “gausswin” function’s
%   default ? of 2.5 and normalized such that the filter integrates to 1,
%   to smooth out any high-frequency changes. 
%
%   The angles are then sequentially checked from head to tail. Every time
%   the angle changes sign or hits 0°, the end of a bend has been found and
%   the count is incremented. Bends found at the start and end of the worm
%   must reflect a segment at least 1/12 the skeleton length in order to be
%   counted. This ignores small bends at the tip of the head and tail.




%JAH NOTE: I have not yet looked at this function ...



% Determine the bend segment length threshold.
lengthThr    = round(size(worm_angles, 1) / 12);
endLengthThr = lengthThr;

% Determine the bend amplitude threshold.
ampThr = 0;

% Are we smoothing the worm angles?
issmoothing = true;


% Should we show the results on the worm(s)?
wormSkeletons = [];
if length(varargin) > 3
    wormSkeletons = squeeze(varargin{4});
end

% Orient the worm correctly.
worm_angles = squeeze(worm_angles);
if size(worm_angles, 1) == 1
    worm_angles = worm_angles';
end

% Compute a guassian filter for the angles.
if issmoothing
    halfLengthThr = round(lengthThr / 2);
    gaussFilter   = gausswin(halfLengthThr * 2 + 1) / halfLengthThr;
end

% Compute the kinks for the worms.
numWorms    = size(worm_angles, 2);
numKinks    = nan(1, numWorms);
for i = 1:size(worm_angles, 2)
    
    % Do we have a worm?
    bends = worm_angles(:,i);
    nanBends = isnan(bends);
    if all(nanBends)
        continue;
    end
    
    % Filter the worm bends.
    if issmoothing
        bends = conv(bends, gaussFilter, 'same');
    end
    
    % Compute the kinks.
    bendSigns = sign(bends);
    kinkAngle = nan(1, length(bendSigns));
    kinkIndex = nan(1, length(bendSigns));
    numKink = 0;
    numSameSign = 0; % the number of adjacent bends with the same sign
    for j = 1:(length(bendSigns) - 1) % the last bend is NaN
        
        % Compute the kink information.
        % Note: data at the zero crossing is counted for both bend sides.
        if bendSigns(j) ~= 0 && ~isnan(bendSigns(j)) && ...
                ~isnan(bendSigns(j + 1)) && ...
                bendSigns(j) ~= bendSigns(j + 1)
            if bendSigns(j) > 0
                amp = max(bends((j - numSameSign):j));
            elseif bendSigns(j) < 0
                amp = min(bends((j - numSameSign):j));
            else
                amp = 0;
            end
            
            % Compute the bend length.
            if numKink == 0
                bendLength = j;
            else
                bendLength = numSameSign + 1;
            end
            
            % Include the zero bend.
            if bendSigns(j+1) == 0
                bendLength = bendLength + 1;
            end
            
            % Add the first kink.
            if numKink == 0
                if bendLength >= endLengthThr && abs(amp) >= ampThr
                    numKink = numKink + 1;
                    kinkAngle(numKink) = amp;
                    kinkIndex(numKink) = j - (bendLength - 1) / 2;
                end
                
                % Add the kink.
            else
                if bendLength >= lengthThr && abs(amp) >= ampThr
                    numKink = numKink + 1;
                    kinkAngle(numKink) = amp;
                    kinkIndex(numKink) = j - (bendLength - 1) / 2;
                end
            end
            
            % Reset the count for adjacent bends with the same sign.
            numSameSign = 0;
            
            % Advance.
        else
            numSameSign = numSameSign + 1;
        end
    end
    
    % Compute the kink information for the last bend.
    bendLength = numSameSign + 1;
    segSigns   = bendSigns((end - bendLength + 1):end);
    segSigns   = segSigns(~isnan(segSigns));
    if isempty(segSigns) || segSigns(1) == 0
        amp = 0;
    elseif segSigns(1) == 1
        amp = max(bends((end - bendLength + 1):end));
    else % if segSigns(1) == -1
        amp = min(bends((end - bendLength + 1):end));
    end
    
    if bendLength >= endLengthThr && abs(amp) >= ampThr
        numKink = numKink + 1;
    end
    
    
    % Record the kink information.
    numKinks(i)    = numKink;
    
    
end
end
