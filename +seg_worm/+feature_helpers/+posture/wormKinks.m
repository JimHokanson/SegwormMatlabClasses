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
%   worm_angles : [48 x n_frames] the worm(s) bend angles at each skeleton point
%
%   Outputs:
%   =======================================================================
%   n_kinks_all : [1 x n_frames] the number of kinks
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
%   default alpha of 2.5 and normalized such that the filter integrates to
%   1, to smooth out any high-frequency changes.
%
%   JAH: Why normalize if no threshold on the amplitude?
%
%   The angles are then sequentially checked from head to tail. Every time
%   the angle changes sign or hits 0°, the end of a bend has been found and
%   the count is incremented. Bends found at the start and end of the worm
%   must reflect a segment at least 1/12 the skeleton length in order to be
%   counted. This ignores small bends at the tip of the head and tail.


% Determine the bend segment length threshold.
length_thr = round(size(worm_angles, 1) / 12); %Threshold is inclusive

% Compute a guassian filter for the angles.
%--------------------------------------------------------------------------
%JAH NOTE: This is a nice way of getting the appropriate odd value
%unlike the other code with so many if statements ...
%- see window code which tries to get an odd value ...
half_length_thr = round(length_thr / 2);
gauss_filter    = gausswin(half_length_thr * 2 + 1) / half_length_thr;


% Compute the kinks for the worms.
n_frames    = size(worm_angles, 2);
n_kinks_all = nan(1, n_frames);
for iFrame = find(any(worm_angles,1))
    
    %Smooth bend angles
    bends = conv(worm_angles(:,iFrame), gauss_filter, 'same');

    n_kinks_all(iFrame) = h__computeNumberOfKinks_New(bends,length_thr);
end


end

function  n_kinks = h__computeNumberOfKinks_New(bends,lengthThr)
    %
    %   We look for sign changes. We get the length of each signed
    %   region and if it is greater than sum length cutoff, we count it as
    %   a kink.
    %
   
    %This code is nearly identical in getForaging
    %-------------------------------------------------------
    n_frames = length(bends);

    dataSign      = sign(bends);
    
    if any(dataSign == 0)
        error('Unhandled code case')
    end
    
    sign_change_I = find(dataSign(2:end) ~= dataSign(1:end-1));

    end_I   = [sign_change_I; n_frames];
    start_I = [1; sign_change_I+1];

    %All NaN values are considered sign changes, remove these ...
    mask = isnan(bends(start_I));
    start_I(mask) = [];
    end_I(mask)   = [];
    
    %The old code had a provision for having NaN values in the middle
    %of the worm. I have not translated that feature to the newer code. I
    %don't think it will ever happen though for a valid frame, only on the
    %edges should you have NaN values.
    if ~isempty(start_I) && any(isnan(bends(start_I(1):end_I(end))))
       error('Unhandled code case')
    end
    %-------------------------------------------------------
    %End of identical code ...
    
    
    lengths = end_I - start_I + 1;
    
    %Adjust lengths for first and last:
    %Basically we allow NaN values to count towards the length for the
    %first and last stretches
    if ~isempty(lengths)
       if start_I(1) ~= 1 %Due to leading NaNs
          lengths(1) = lengths(1) + start_I(1)-1;  
       end
       if end_I(end) ~= n_frames %Due to trailing NaNs
          lengths(end) = lengths(end) + (n_frames - end_I(end));
       end
    end
    
    n_kinks = sum(lengths >= lengthThr);
   
end
