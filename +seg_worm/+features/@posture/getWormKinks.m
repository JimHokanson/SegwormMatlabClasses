function getWormKinks(obj, bend_angles,p_opts)
%
%
%   seg_worm.features.posture.getWormKinks
%
%   Old Name: 
%   - wormKinks.m
%
%   Inputs:
%   =======================================================================
%   bend_angles : [48 x n_frames] the worm(s) bend angles at each skeleton point
%
%   Populates:
%   =======================================================================
%   kinks
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
%   JAH: No normalization needed as there is currently no amplitude
%   threshold and we're only counting #, not size
%
%   The angles are then sequentially checked from head to tail. Every time
%   the angle changes sign or hits 0°, the end of a bend has been found and
%   the count is incremented. Bends found at the start and end of the worm
%   must reflect a segment at least 1/12 the skeleton length in order to be
%   counted. This ignores small bends at the tip of the head and tail.

%{

Missing Old Features
-------------------------------------
1) Thresholding the amplitude of a bend to determine whether or not to
count it. No amplitude threshold was in place so this code was removed.
Currently any sign change is considered a bend as long as it is long
enough.


%}


% Determine the bend segment length threshold.
n_angles = size(bend_angles, 1);
length_threshold = round(n_angles * p_opts.KINK_LENGTH_THRESHOLD_PCT);

% Compute a gaussian filter for the angles.
%--------------------------------------------------------------------------
%JAH NOTE: This is a nice way of getting the appropriate odd value
%unlike the other code with so many if statements ...
%- see window code which tries to get an odd value ...
%- I'd like to go back and fix that code ...
half_length_thr = round(length_threshold / 2);
gauss_filter    = gausswin(half_length_thr * 2 + 1) / half_length_thr;


% Compute the kinks for the worms.
n_frames    = size(bend_angles, 2);
n_kinks_all = NaN(1, n_frames);
for iFrame = find(any(bend_angles,1))
    
    smoothed_bend_angles = conv(bend_angles(:,iFrame), gauss_filter, 'same');

    n_kinks_all(iFrame) = h__computeNumberOfKinks_New(smoothed_bend_angles,length_threshold);
end

obj.kinks = n_kinks_all;


end

% % function h__explainFunction(orig_worm_bend_angles,smoothed_bend_angles,sign_change_I,n_bends)
% %    %Not sure how I want to get access to this information ...
% %    %
% %    %Would also be nice to have access to the original worm for plotting the
% %    %skeleton as well ...
% %    %
% %    %Would like to avoid code duplication ...
% % end

function  n_kinks = h__computeNumberOfKinks_New(smoothed_bend_angles,length_threshold)
    %
    %   We look for sign changes. We get the length of each signed
    %   region and if it is greater than sum length cutoff, we count it as
    %   a kink.
    %
    %   Inputs
    %   ===================================================================
    %   smoothed_bend_angles : [1 x 48]
    %   length_threshold     : scalar, 
   
    %This code is nearly identical in getForaging
    %-------------------------------------------------------
    n_frames = length(smoothed_bend_angles);

    dataSign      = sign(smoothed_bend_angles);
    
    if any(dataSign == 0)
        %I don't expect that we'll ever actually reach 0
        %The code for zero was a bit weird, it keeps counting if no sign
        %change i.e. + + + 0 + + + => all +
        %
        %but if counts for both if sign change
        % + + 0 - - - => 3 +s and 4 -s
        error('Unhandled code case')
    end
    
    sign_change_I = find(dataSign(2:end) ~= dataSign(1:end-1));

    end_I   = [sign_change_I; n_frames];
    start_I = [1; sign_change_I+1];

    %All NaN values are considered sign changes, remove these ...
    mask = isnan(smoothed_bend_angles(start_I));
    start_I(mask) = [];
    end_I(mask)   = [];
    
    %The old code had a provision for having NaN values in the middle
    %of the worm. I have not translated that feature to the newer code. I
    %don't think it will ever happen though for a valid frame, only on the
    %edges should you have NaN values.
    if ~isempty(start_I) && any(isnan(smoothed_bend_angles(start_I(1):end_I(end))))
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
    
    n_kinks = sum(lengths >= length_threshold);
       
end
