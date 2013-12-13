function n_kinks_all = wormKinks(worm_angles)
%
%
%   seg_worm.feature_helpers.posture.wormKinks
%
%   Old Name: wormKinks.m
%
%
%   Inputs:
%   =======================================================================
%   worm_angles : the worm(s) bend angles at each skeleton point
%
%   Outputs:
%       numKinks - the number of kinks
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
%   JAH: Why normalize if no threshold ?????
%
%   The angles are then sequentially checked from head to tail. Every time
%   the angle changes sign or hits 0°, the end of a bend has been found and
%   the count is incremented. Bends found at the start and end of the worm
%   must reflect a segment at least 1/12 the skeleton length in order to be
%   counted. This ignores small bends at the tip of the head and tail.


% Determine the bend segment length threshold.
lengthThr = round(size(worm_angles, 1) / 12); %Inclusive threshold

% Compute a guassian filter for the angles.
%--------------------------------------------------------------------------
%JAH NOTE: This is a nice way of getting the appropriate odd value
%unlike the other code with so many if statements ...
%- see window code which tries to get an odd value ...
halfLengthThr = round(lengthThr / 2);
gaussFilter   = gausswin(halfLengthThr * 2 + 1) / halfLengthThr;


% Compute the kinks for the worms.
n_frames    = size(worm_angles, 2);
n_kinks_all    = nan(1, n_frames);
for iFrame = find(any(worm_angles,1))
    
    bends = conv(worm_angles(:,iFrame), gaussFilter, 'same');

    n_kinks_all(iFrame) = h__computeNumberOfKinks_New(bends,lengthThr);
end


end

function  n_kinks = h__computeNumberOfKinks_New(bends,lengthThr)
    %
    %   We look for sign changes. We get the length of each signed
    %   region and if it is greater than sum length cutoff, we count it as
    %   a kink.
    %
   
   %TODO: Check for 0 sign, this is unhandled ...
   %TODO: Should check for any NaN between start_I(1) and end_I(end)
   % - this should never happpen
   
    %This code is identical in getForaging
    %-------------------------------------------------------
    n_frames = length(bends);

    dataSign      = sign(bends);
    sign_change_I = find(dataSign(2:end) ~= dataSign(1:end-1));

    end_I   = [sign_change_I; n_frames];
    start_I = [1; sign_change_I+1];

    %All NaN values are considered sign changes, remove these ...
    mask = isnan(bends(start_I));
    start_I(mask) = [];
    end_I(mask)   = [];
    %-------------------------------------------------------
    %End of identical code ...
    
    
    lengths = end_I - start_I + 1;
    
    %Adjust lengths for first and last:
    %Basically we allow NaN values to count towards the length for the
    %first and last
    if ~isempty(lengths)
       if start_I(1) ~= 1 %Due to NaN
          lengths(1) = lengths(1) + start_I(1)-1;  
       end
       if end_I(end) ~= n_frames
          lengths(end) = lengths(end) + (n_frames - end_I(end));
       end
    end
    
    n_kinks = sum(lengths >= lengthThr);
   
end
